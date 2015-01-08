
from rootpy.plotting import Hist
from rootpy.plotting.views import ScaleView
from rootpy.io import File
try:
    from collections import OrderedDict
except ImportError:
    from ordered import OrderedDict

##helper classes
def xfrange(start, stop, step):
    while start < stop:
        yield start
        start += step

def getColorList(ncolors):
    from colorsys import hls_to_rgb,rgb_to_hls
    from math import pi

    from random import random
    colors=[]
    for i in xfrange(0,2.*pi,(2.*pi)/ncolors):
        #colorsys.hls_to_rgb(h, l, s)
        #colors.append(hls_to_rgb(i, random() , pi/2. + random() ))
        colors.append(hls_to_rgb(i, random() , random() ))
    return colors

def duke_errorbar(hists,
             xerr=True, yerr=True,
             xpadding=0, ypadding=.1,
             xerror_in_padding=True,
             yerror_in_padding=True,
             emptybins=True,
             snap=True,
             axes=None,
             ignore_binns=None,
             **kwargs):
    from rootpy.plotting.hist import _Hist
    from rootpy.plotting.graph import _Graph1DBase
    from rootpy.plotting.root2matplotlib import _set_bounds

    """
    Make a matplotlib errorbar plot from a ROOT histogram or graph
    or list of histograms and graphs.

    Parameters
    ----------

    hists : Hist, Graph or list of Hist and Graph
        The histogram(s) and/or Graph(s) to be plotted

    xerr : bool, optional (default=True)
        If True, x error bars will be displayed.

    yerr : bool or string, optional (default=True)
        If False, no y errors are displayed.  If True, an individual y
        error will be displayed for each hist in the stack.  If 'linear' or
        'quadratic', a single error bar will be displayed with either the
        linear or quadratic sum of the individual errors.

    xpadding : float or 2-tuple of floats, optional (default=0)
        Padding to add on the left and right sides of the plot as a fraction of
        the axes width after the padding has been added. Specify unique left
        and right padding with a 2-tuple.

    ypadding : float or 2-tuple of floats, optional (default=.1)
        Padding to add on the top and bottom of the plot as a fraction of
        the axes height after the padding has been added. Specify unique top
        and bottom padding with a 2-tuple.

    xerror_in_padding : bool, optional (default=True)
        If True then make the padding inclusive of the x errors otherwise
        only pad around the x values.

    yerror_in_padding : bool, optional (default=True)
        If True then make the padding inclusive of the y errors otherwise
        only pad around the y values.

    emptybins : bool, optional (default=True)
        If True (the default) then plot bins with zero content otherwise only
        show bins with nonzero content.

    snap : bool, optional (default=True)
        If True (the default) then the origin is an implicit lower bound of the
        histogram unless the histogram has both positive and negative bins.

    axes : matplotlib Axes instance, optional (default=None)
        The axes to plot on. If None then use the global current axes.

    kwargs : additional keyword arguments, optional
        All additional keyword arguments are passed to matplotlib's errorbar
        function.

    Returns
    -------

    The return value from matplotlib's errorbar function, or list of such
    return values if a list of histograms and/or graphs was plotted.

    """
    if axes is None:
        axes = plt.gca()
    curr_xlim = axes.get_xlim()
    curr_ylim = axes.get_ylim()
    was_empty = not axes.has_data()
    if isinstance(hists, (_Hist, _Graph1DBase)):
        # This is a single plottable object.
        returns = _duke__errorbar(
            hists, xerr, yerr,
            axes=axes, emptybins=emptybins,
            ignore_binns=ignore_binns,
             **kwargs)
        _set_bounds(hists, axes=axes,
                    was_empty=was_empty,
                    prev_ylim=curr_ylim,
                    xpadding=xpadding, ypadding=ypadding,
                    xerror_in_padding=xerror_in_padding,
                    yerror_in_padding=yerror_in_padding,
                    snap=snap)
    else:
        returns = []
        for h in hists:
            returns.append(duke_errorbar(
                h, xerr=xerr, yerr=yerr, axes=axes,
                xpadding=xpadding, ypadding=ypadding,
                xerror_in_padding=xerror_in_padding,
                yerror_in_padding=yerror_in_padding,
                snap=snap,
                emptybins=emptybins,
                **kwargs))
    return returns


def _duke__errorbar(h, xerr, yerr, axes=None, emptybins=True, ignore_binns=None, zorder=None, **kwargs):
    from rootpy.plotting.root2matplotlib import _set_defaults
    import numpy as np
    if axes is None:
        axes = plt.gca()
    if zorder is None:
        zorder = _get_highest_zorder(axes) + 1
    _set_defaults(h, kwargs, ['common', 'errors', 'errorbar', 'marker'])
    if xerr:
        xerr = np.array([list(h.xerrl()), list(h.xerrh())])
    if yerr:
        yerr = np.array([list(h.yerrl()), list(h.yerrh())])
    x = np.array(list(h.x()))
    y = np.array(list(h.y()))
    remove=[]
    if ignore_binns is not None:
        for i in range(len(x)-1,-1,-1):
            ignore=False
            for ihist in ignore_binns:
                if ihist[i+1].value==0:
                    ignore=True
            if ignore:
                remove.append(i)
    x=np.delete(x, remove)
    y=np.delete(y, remove)
    if xerr is not False:
        xerr=np.delete(xerr, remove,1)
    if yerr is not False:
        yerr=np.delete(yerr, remove,1)

    if not emptybins:
        nonempty = y != 0
        x = x[nonempty]
        y = y[nonempty]
        if xerr is not False:
            xerr = xerr[:, nonempty]
        if yerr is not False:
            yerr = yerr[:, nonempty]
    return axes.errorbar(x, y, xerr=xerr, yerr=yerr, zorder=zorder, **kwargs)





##@class HistSorage
# Class to hangle histograms functions
#
#
#
# written by Klaas Padeken 2015
class HistSorage(object):
    def __init__(self,xs, lumi, path=None,isData=False):
        self.views=OrderedDict()
        self.hists=OrderedDict()
        self.files=OrderedDict()
        self.genNumber={}
        self.basepath=path
        self.verbosity=3
        self._scaled=False
        self.datadrivenHist=""
        self.xs=xs
        self.lumi=lumi
        self.isData=isData
        self._joinList=False

    def __del__(self):
        for name in  self.files:
            self.files[name].Close()

    def _getGenNumbers(self):
        for name in self.files:
            if name in self.genNumber or self.datadrivenHist==name:
                continue
            try:
                counter=self.files[name].Get("h_counters")
                Nev=max(counter[1].value,counter[2].value)
            except:
                if self.verbosity==3:
                    print("[Info] If you want to use Nev from the root file store 'h_counters'")
                    print("will set to 1 for %s"%(name))
                Nev=1
            self.genNumber[name]=Nev

    def _addToScaledView(self):
        for name in self.files:
            if name in self.views:
                continue
            if name==self.datadrivenHist or self.isData:
                weight=1.
            else:
                weight=self.xs[name].as_float("xs")*self.xs[name].as_float("weight")*self.lumi/self.genNumber[name]
            self.views[name]=ScaleView(self.files[name],weight)
            self._scaled=True

    def addAllFiles(self,tag="",veto=None):
        if self.basepath==None:
            raise RuntimeError("You must set a basepath to add all files from one directory!")
        import glob
        fileList=glob.glob(self.basepath+"/*"+tag+"*.root")
        for file in fileList:
            if veto is not None:
                vetoed=False
                for v in veto:
                    if v.lower() in file.split("/")[-1].lower():
                        vetoed=True
                if vetoed:
                    continue
            name=file.split("/")[-1].replace(".root","")
            self.files[name]=File(file, "read")
        self._getGenNumbers()
        self._addToScaledView()

    def addFile(self,name):
        self.files[name]=File(self.basepath+"/"+name+".root", "read")
        self._getGenNumbers()
        self._addToScaledView()

    def addFileList(self,fileList):
        if type(fileList)==type(list()):
            for file in fileList:
                self.files[file]=File(self.basepath+"/"+file+".root", "read")
        if isinstance(fileList,dict):
            import itertools
            useList=list(itertools.chain.from_iterable(fileList.values()))
            for file in useList:
                self.files[file]=File(self.basepath+"/"+file+".root", "read")
            self._joinList=fileList
        self._getGenNumbers()
        self._addToScaledView()


    def addJoinList(self,joinList):
        self._joinList=joinList

    def clearHists(self):
        self.hists={}

    def getAdded(self,name):
        temp = []
        for key in self.hists:
            if name in key:
                temp.append( self.hists[key] )
        return sum(temp)

    def getAllAdded(self,ignoreScale=False):
        if not self._scaled and not ignoreScale:
            raise RuntimeError("Add all histograms without scaling. I think not!")
        return sum(self.hists.values())

    def getHist(self,hist):
        for f in self.views:
            self.hists[f]=self.views[f].Get(hist)
            self.hists[f].Sumw2()
        if self._joinList is not False:
            self.joinList(self._joinList)

    def join(self,name,label):
        if name == "*" or name == "":
            # join all bg:
            joined = self.getAllAdded()
            self.hists=OrderedDict()
            self.hists[label]=joined
        else:
            # join hist filtered by name:
            joined = self.getAdded(name)
            for key in self.hists:
                if name in key:
                    self.hists.pop(key)
            self.hists[label]=joined

    def joinList(self,joinList):
        for name in joinList:
            self.hists[name]=self.hists[joinList[name][0]]
            self.hists.pop(joinList[name][0])
            for h in joinList[name][1:]:
                self.hists[name]+=self.hists[h]
                self.hists.pop(h)

    def rebin(self,width=0,factor=0):
        if width!=0:
            factor=int(width/self.hists.values()[-1].xwidth(1)+0.5)
        for name in self.hists:
            self.hists[name].Rebin(factor)

    def scale_cfg(self):
        if self._scaled==True:
            raise RuntimeError("Histograms already scaled!")
        for name in self.hists:
            if name==self.datadrivenHist:
                continue
            weight=self.xs[name].as_float("xs")*self.xs[name].as_float("weight")*self.lumi/self.genNumber[name]
            self.hists[name].Scale(weight)
        self._scaled=True

    def setStyle(self,style="bg",colors=None):
        if style=="bg":
            for key in self.hists:
                self.hists[key].fillstyle = 'solid'
                self.hists[key].linewidth = 0
        if style=="sg":
            for key in self.hists:
                self.hists[key].fillstyle = '0'
                self.hists[key].linewidth = 1
        for key in self.hists:
            self.hists[key].xaxis.SetTitle("$\mathsf{"+self.hists[key].xaxis.GetTitle()+"}$")
            self.hists[key].yaxis.SetTitle("Events/%.f GeV"%(self.hists.values()[-1].xwidth(1)))
            self.hists[key].SetTitle(key)
            if colors!=None:
                if isinstance(colors, (list)):
                    usecolor=colors.pop()
                    self.hists[key].fillcolor = usecolor
                    self.hists[key].linecolor = usecolor
                if isinstance(colors, (dict)):
                    self.hists[key].fillcolor = colors[key]
                    self.hists[key].linecolor = colors[key]

    def setDataDriven(self,ddhist):
        self.datadrivenHist=ddhist

    def getHistList(self):
        return self.hists.values()
















#class Hist(object):
    #def __init__(self,h,name=""):
        #if name=="":
            #self.name=h.GetName()
        #else:
            #self.name=name
        #self.hist=h
        #self.fcolor=1
        #self.lcolor=1
        #self.lstyle=1
        #self.lwidth=3
        #self.fstyle=0
        #self.label=""
        #self.noLine=False

        #self.scale=1.
        #self._scaled=False
        #self._EventNumber=0
        #uniqueName(h)

    #def update(self):
        #self.hist.UseCurrentStyle()
        ##lines
        #self.hist.SetLineStyle(self.lstyle)
        #self.hist.SetLineColor(self.lcolor)
        #self.hist.SetLineWidth(self.lwidth)

        ##fill
        #self.hist.SetFillStyle(self.fstyle)
        #self.hist.SetFillColor(self.fcolor)
        #if self.noLine:
            #self.hist.SetLineColor(self.fcolor)

    #def scale(self,force=False,sf=None):
        #if sf is None:
            #sf=self.scale
        #if self._scaled and not forced:
            #raise RuntimeError("You can not scale a already scaled hist!\n Use forced=True if you must!!")
        #self.hist.Scale(sf)
        #self._scaled=True

    #def Scale(self,sf):
        #scale(sf=sf)

    #def setEventNumber(self,n):
        #self._EventNumber=n

    #def add(self,other):
        #h=Hist(self.hist,name=self.name+"_"+other.name)
        #h.hist.Add(other.hist)
        #return h

    #def __add__(self,other):
        #return self.add(other)

    #def Add(self,other):
        #return self.add(other)

    #def __str__(self):
        #return str(self.fcolor)+" "+str(self.lcolor)+" "+str(self.lstyle)+" "+str(self.lwidth)+" "+str(self.fstyle)+" "+str(self.scale)



#class Plot(object):
    #def __init__(self,name):
        ##public
        #self.bg, self.sg = [], []
        #self.bgDict,self.sgDict = {}, {}
        #self.data=None
        #self.name=name
        #self.fileName=""
        #self.rangex=[]
        #self.rangey=[]
        #self.rebin=1
        #self.variable=False
        #self.binning=[]
        #self.useFixedWidth=0
        #self.smallestBinWidth=-1
        #self.yAxisLabel="Events / %s "
        #self.yAxisUnit=""
        ##self.legend=Legend()

        #self.canvas=canvas=ROOT.TCanvas("c","c",1050,1050)
        #self.legend=rootplotlib.Legend(pad=ROOT.gPad)
        #rootplotlib.init()
        #self.canvas.UseCurrentStyle()

        ##private
        #self._knownUnits=["GeV","TeV","rad"]

    #def draw(self,pad=None):
        #"""draw the standard stacked hist plot"""
        #if pad is None:
            #pad=self.canvas

        #first=None
        ##make bg if there

        #hs=None
        #if len(self.bg)>0:
            #hs=ROOT.THStack("hs","hs")
            #for b in self.bg:
                #hs.Add(b.hist)
            #if first is None:
                #first=hs
                #hs.Draw("hist")
            #else:
                #hs.Draw("hist same")
            #hs.GetXaxis().SetTitle(self.bg[-1].hist.GetXaxis().GetTitle())
            #hs.GetYaxis().SetTitle(self.bg[-1].hist.GetYaxis().GetTitle())

        #if len(self.sg)>0:
            #for s in self.sg:
                #if first is None:
                    #first=s
                    #s.hist.Draw("hist")
                #else:
                    #s.hist.Draw("same hist")

        #if self.data is not None:
            #if first is None:
                #first=self.data.hist
                #self.data.hist.Draw("EP")
            #else:
                #self.data.hist.Draw("EP same")
        #ROOT.gPad.RedrawAxis("g")
        #ROOT.gPad.Update()
        #pad.SaveAs("test.png")
        #raw_input("Nnk")


    #def applyStyle(self):
        #"""apply all the plot styles one can set on every hist"""

        ##fist the hists itself
        #all_histsW=[i for i in (self.bg+self.sg)]
        #all_histsW+=[self.data]

        #for h in all_histsW:
            #h.update()

        ##now the plot
        #all_hists=[i.hist for i in all_histsW]
        #for h in all_hists:
            #if len(self.rangex)==2:
                #h.GetXaxis().SetRangeUser(self.rangex[0],self.rangex[1])
            #if len(self.rangey)==2:
                #h.GetYaxis().SetRangeUser(self.rangey[0],self.rangey[1])
            #if self.variable:
                #h=self.variableBinning(h)
            #else:
                #h.Rebin(self.rebin)
            #self.setYaxisLabel(h)

        #for b in self.bg:
            #self.legend.add( b.hist, b.label, "f")
        #for s in self.sg:
            #self.legend.add( s.hist, s.label, "l")



    #def variableBinning(self, h):
        #if self.variable==False or len(self.binning)==0:
            #raise ValueError('varible binning is not se')

        #smallestBinW=-1
        #for i in self.binning:
            #if self.useFixedWidth:
                #smallestBinW = self.smallestBinWidth
            #elif i[0]<smallestBinW or smallestBinW<0:
                #smallestBinW=i[0]
            #for ibin in xfrange(i[1],i[2],i[0]):
                #binVector.append(ibin)
        #binVector.append(binVector[-1]+binning[-1][0])
        #htmp=h.Rebin(len(binVector)-1,h.GetName()+str(id(h)),array("d",binVector))

        #self.smallestBinWidth=smallestBinW

        #for i in range(htmp.GetNbinsX()+1):
            #sf = smallestBinW/(htmp.GetBinWidth(i))
            #htmp.SetBinContent(i, sf*(htmp.GetBinContent(i)))
            #htmp.SetBinError(i, sf*(htmp.GetBinError(i)))
        #return  htmp

    #def setYaxisLabel(self, h):
        #from rounding import rounding
        #if self.variable:
            #binWidth=self.smallestBinWidth
        #else:
            ##wtf why does this not work root version??
            #binWidth=h.GetXaxis().GetBinWidth(1)
        #roundObj=rounding()
        #yax=self.yAxisLabel%(roundObj.latex(binWidth))
        #if self.yAxisUnit!="":
            #yax+=self.yAxisUnit
        #else:
            #xax=h.GetXaxis().GetTitle()
            #for unit in self._knownUnits:
                #if unit in xax:
                    #yax+=unit
        #h.GetYaxis().SetTitle(yax)

    #def addBG(self,h,name):
        #h.fstyle=1001
        #h.noLine=True

        #self.bgDict.update({name:h})
        #self.bg.append(h)

    #def addSG(self,h,name):
        #self.sgDict.update({name:h})
        #self.sg.append(h)

    #def addData(self,h):
        #self.data=h

    #def getAllBGAdded(self):
        #addedBG=self.bg[0]
        #for i in self.bg[1:]:
            #addedBG += i
        #addedBG.name="allbg"
        #return addedBG

    #def getBG(self,name):
        #for b in self.bg:
            #if b.name==name:
                #return b
        #return None

    #def getBGAdded(self,name,label):
        #temp = []
        #for key in self.bgDict.keys():
            #if name in key:
                #temp.append( self.bgDict[key] )
        #addedBG = temp[0]
        #for i in temp[1:]:
            #addedBG += i
        #addedBG.name=label
        #return addedBG

    #def getSG(self,name):
        #for s in self.sg:
            #if s.name==name:
                #return s
        #print "SG not found"
        #return None

    #def joinBG(self,name,label):
        #if name == "*" or name == "":
            ## join all bg:
            #joined = self.getAllBGAdded()
            #self.bg = []
            #self.bgDict = {}
            #self.addBG(joined, label)
        #else:
            ## join bg filtered by name:
            #joined = self.getBGAdded(name,label)

            ## delete bg from plot dictionary (todo: also delete from bg list):
            #for key in self.bgDict.keys():
                #if name in key:
                    #self.bgDict.pop(key,0)

            #self.addBG(joined, label)


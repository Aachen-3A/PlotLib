

from rootpy.plotting import Hist
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




##@class HistSorage
# Class to hangle histograms functions
#
#
#
# written by Klaas Padeken 2015
class HistSorage(object):
    def __init__(self,path=None):
        self.files=OrderedDict()
        self.hists=OrderedDict()
        self.genNumber={}
        self.basepath=path
        self.verbosity=3
        self._scaled=False

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

    def addFile(self,name):
        self.files[name]=File(self.basepath+"/"+name+".root", "read")

    def getHist(self,hist):
        for f in self.files:
            self.hists[f]=self.files[f].Get(hist)
            self.hists[f].Sumw2()
        if self.genNumber=={}:
            self._getGenNumbers()

    def scale_cfg(self,lumi,xs):
        #from configobj import ConfigObj
        for name in self.hists:
            weight=xs[name].as_float("xs")*xs[name].as_float("weight")*lumi/self.genNumber[name]
            self.hists[name].Scale(weight)
        self._scaled=True

    def rebin(self,width=0,factor=0):
        if width!=0:
            factor=int(width/self.hists.values()[-1].xwidth(1)+0.5)
        for name in self.hists:
            self.hists[name].Rebin(factor)

    def clearHists(self):
        self.hists={}

    def getAllAdded(self,ignoreScale=False):
        if not self._scaled and not ignoreScale:
            raise RuntimeError("Add all histograms without scaling. I think not!")
        return sum(self.hists.values())

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

    def getAdded(self,name):
        temp = []
        for key in self.hists:
            if name in key:
                temp.append( self.hists[key] )
        return sum(temp)

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
            self.hists[key].xaxis.SetTitle("$\mathrm{"+self.hists[key].xaxis.GetTitle()+"}$")
            self.hists[key].yaxis.SetTitle("Events/%.f GeV"%(self.hists.values()[-1].xwidth(1)))
            self.hists[key].SetTitle(key)
            if colors!=None:
                usecolor=colors.pop()
                self.hists[key].fillcolor = usecolor
                self.hists[key].linecolor = usecolor

    def getHistList(self):
        return self.hists.values()

    def _getGenNumbers(self):
        for name in self.files:
            try:
                counter=self.files[name].Get("h_counters")
                Nev=max(counter[1].value,counter[2].value)
            except:
                if self.verbosity==3:
                    print("[Info] If you want to use Nev from the root file store 'h_counters'")
                    print("will set to 1")
                Nev=1
            self.genNumber[name]=Nev














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


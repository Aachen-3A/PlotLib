
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

def getDictValue(hist,parmDict):
    for par in parmDict:
        if par in hist:
            return parmDict[par]
    return None


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
# To store a list of hists, which are scaled, joined, rebinned and otherwise
# manipulated.
#
# written by Klaas Padeken 2015
class HistSorage(object):
    ## Init function
    #
    # In this function the default variables are set and initialized.
    # @param[in] xs is a xs ConfigObj needed for scaling
    # @param[in] lumi is the lumi in pb
    # @param[in] path is the default path of the files (default=None)
    # @param[in] isData is a switch  (default=None)
    def __init__(self,xs, lumi, path=None,isData=False):
        self.views=OrderedDict()
        self.hists=OrderedDict()
        self.files=OrderedDict()
        self.colorList={}
        self.genNumber={}
        self.basepath=path
        self.verbosity=3
        self._scaled=False
        self.datadrivenHist=""
        self.xs=xs
        self.lumi=lumi
        self.isData=isData
        self._joinList=False

    ## del function
    #
    # This deletes the main objects nedded to not get a crash at the end!
    def __del__(self):
        for name in  self.files:
            self.files[name].Close()
    ##------------------------------------------------------------------
    ## Private functions
    ##------------------------------------------------------------------
    ## Function to get the event numbers from h_counter in file!
    #
    # The function fills the dict genNumber with the event numbers.
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

    ## Function to add files to a scaled view.
    #
    # The scaled view dict "views" now retruns all histograms scaled!
    def _addToScaledView(self):
        for name in self.files:
            if name in self.views:
                continue
            if name==self.datadrivenHist or self.isData:
                weight=1.
            else:
                if "weight" in self.xs[name]:
                    weight=self.xs[name].as_float("xs")*self.xs[name].as_float("weight")*self.lumi/self.genNumber[name]
                else:
                    weight=self.xs[name].as_float("xs")*self.lumi/self.genNumber[name]
            self.views[name]=ScaleView(self.files[name],weight)
            self._scaled=True

    ## Function to add all files in the given path
    #
    # Use setPath(path) to set the path if did not in the init.
    # @param[in] tag if regexpr is not used all *.root files containing the tag are added
    # @param[in] veto define a !list!! of veto strings not case sensitive
    # @param[in] regexpr use a regular expression to find the file names (need .root at the end if
    # you want to use root files!!
    # @param[in] joinName if specified all files matching the expressions above will be added to the list of files that should be joined.
    def addAllFiles(self,tag="",veto=None, regexpr=None, joinName=None):
        if self.basepath==None:
            raise RuntimeError("You must set a basepath to add all files from one directory!")
        if regexpr is not None:
            import re,os
            fileList = [f for f in os.listdir(self.basepath+"/") if re.search(r'%s'%(regexpr), f)]
        else:
            import glob
            fileList=glob.glob(self.basepath+"/*"+tag+"*.root")
        tmpList=[]
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
            tmpList.append(self.files[name])
        self._getGenNumbers()
        self._addToScaledView()
        if joinName is not None:
            if self._joinList is not False:
                self._joinList[joinName]=tmpList
            else:
                self._joinList=OrderedDict()
                self._joinList[joinName]=tmpList


    ## Function to add a single file
    #
    # Use setPath(path) to set the path if did not in the init.
    # @param[in] name the name of the file that should be added!
    def addFile(self,name):
        self.files[name]=File(self.basepath+"/"+name+".root", "read")
        self._getGenNumbers()
        self._addToScaledView()

    ## Function to add files specified as a list or (ordered)dict
    #
    # Use setPath(path) to set the path if did not in the init.
    # @param[in] fileList list or dict of the files you want to add
    # if the dict is used the files are joined to a single hist with this key
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

    ## Function to add a dict to join files
    #
    # @param[in] joinList wich should be a (ordered)dict
    def addJoinList(self,joinList):
        self._joinList=joinList

    ## Function to clear hists
    #
    # use this if you want to plot a new set of hists
    def clearHists(self):
        self.hists=OrderedDict()

    ## Function to get a hist that is the sum of hists in the storage
    #
    # handy if you want only a subgroup as a hist
    # @param[in] name add only files that contain the name (default="")
    # @param[in] ignoreScale if you want to add hists that are not scaled (default=False)
    # @param[out] Hist
    def getAdded(self,name="",ignoreScale=False):
        if not self._scaled and not ignoreScale:
            raise RuntimeError("Add all histograms without scaling. I think not!")
        temp = []
        for key in self.hists:
            if name in key:
                temp.append( self.hists[key] )
        return sum(temp)

    ## Function to get a hist that is the sum of all hists in the storage
    #
    # same as getAdded() perhaps faster
    # @param[in] ignoreScale if you want to add hists that are not scaled (default=False)
    # @param[out] Hist
    def getAllAdded(self,ignoreScale=False):
        if not self._scaled and not ignoreScale:
            raise RuntimeError("Add all histograms without scaling. I think not!")
        return sum(self.hists.values())

    ## Function get hists from files
    #
    # the hists ate added to .hists and joined if a joinList exist
    # @param[in] ignoreScale if you want to add hists that are not scaled (default=False)
    # @param[out] Hist
    def getHist(self,hist):
        for f in self.views:
            self.hists[f]=self.views[f].Get(hist)
            self.hists[f].Sumw2()
        if self._joinList is not False:
            self.joinList(self._joinList)

    ## Function join files containing name to one nameed label
    #
    # @param[in] name add all files containing name
    # @param[in] label name of the resulting new hist
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

    ## Function join files via a (ordered)dict
    #
    # @param[in] joinList add all files that are in the (ordered)dict to one hist with the name of the key
    def joinList(self,joinList):
        for name in joinList:
            self.hists[name]=self.hists[joinList[name][0]]
            self.hists.pop(joinList[name][0])
            for h in joinList[name][1:]:
                self.hists[name]+=self.hists[h]
                self.hists.pop(h)

    ## Function rebin the all hists
    #
    # @param[in] width try to rebin to a specific width
    # @param[in] factor rebin to with a factor
    # if both are given the width is used
    def rebin(self,width=0,factor=0):
        if width!=0:
            factor=int(width/self.hists.values()[-1].xwidth(1)+0.5)
        for name in self.hists:
            self.hists[name].Rebin(factor)

    ## Function to set the datadiven name flag
    #
    # @param[in] ddhist the name of the datadriven hist
    def setDataDriven(self,ddhist):
        self.datadrivenHist=ddhist

    ## Function to set path of the hists
    #
    # @param[in] path
    def setPath(self,path):
        self.basepath=path

    ## Function to set the style of the histograms
    #
    # @param[in] style "bg" and "sg" posible (default="bg")
    # @param[in] colors a list/dict of colors that the hists should have
    # if colors is not specified the internal colorListis used if set

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
                    if key in colors:
                        self.hists[key].fillcolor = colors[key]
                        self.hists[key].linecolor = colors[key]
            elif len(self.colorList)>0:
                if key in self.colorList:
                    self.hists[key].fillcolor = self.colorList[key]
                    self.hists[key].linecolor = self.colorList[key]

    ## Function to get the hists as a list
    #
    # @param[out] list of all stored hists
    def getHistList(self):
        return self.hists.values()




import ROOT
import rootplotlib


##helper classes
def xfrange(start, stop, step):
    while start < stop:
        yield start
        start += step
def uniqueName(h):
    h.SetName(h.GetName()+str(id(h)))
    
    


class Hist(object):
    def __init__(self,h,name=""):
        if name=="":
            self.name=h.GetName()
        else:
            self.name=name
        self.hist=h
        self.fcolor=1
        self.lcolor=1
        self.lstyle=1
        self.lwidth=3
        self.fstyle=0
        self.label=""
        self.noLine=False
        
        self.scale=1.
        self._scaled=False
        self._EventNumber=0
        uniqueName(h)
    
    def update(self):
        self.hist.UseCurrentStyle()
        #lines
        self.hist.SetLineStyle(self.lstyle)
        self.hist.SetLineColor(self.lcolor)
        self.hist.SetLineWidth(self.lwidth)
        
        #fill
        self.hist.SetFillStyle(self.fstyle)
        self.hist.SetFillColor(self.fcolor)
        if self.noLine:
            self.hist.SetLineColor(self.fcolor)
    
    def scale(self,force=False,sf=None):
        if sf is None:
            sf=self.scale
        if self._scaled and not forced:
            raise RuntimeError("You can not scale a already scaled hist!\n Use forced=True if you must!!")
        self.hist.Scale(sf)
        self._scaled=True
    
    def Scale(self,sf):
        scale(sf=sf)
        
    def setEventNumber(self,n):
        self._EventNumber=n
    
    def add(self,other):
        h=Hist(self.hist,name=self.name+"_"+other.name)
        h.hist.Add(other.hist)
        return h
    
    def __add__(self,other):
        return self.add(other)
    
    def Add(self,other):
        return self.add(other)
    
    def __str__(self):
        return str(self.fcolor)+" "+str(self.lcolor)+" "+str(self.lstyle)+" "+str(self.lwidth)+" "+str(self.fstyle)+" "+str(self.scale)



class Plot(object):
    def __init__(self,name):
        #public
        self.bg, self.sg = [], []
        self.bgDict,self.sgDict = {}, {}
        self.data=None
        self.name=name
        self.fileName=""
        self.rangex=[]
        self.rangey=[]
        self.rebin=1
        self.variable=False
        self.binning=[]
        self.useFixedWidth=0
        self.smallestBinWidth=-1
        self.yAxisLabel="Events / %s "
        self.yAxisUnit=""
        #self.legend=Legend()
        
        self.canvas=canvas=ROOT.TCanvas("c","c",1050,1050)
        self.legend=rootplotlib.Legend(pad=ROOT.gPad)
        rootplotlib.init()
        self.canvas.UseCurrentStyle()

        #private
        self._knownUnits=["GeV","TeV","rad"]
        
    def draw(self,pad=None):
        """draw the standard stacked hist plot"""
        if pad is None:
            pad=self.canvas
        
        first=None
        #make bg if there
        
        hs=None
        if len(self.bg)>0:
            hs=ROOT.THStack("hs","hs") 
            for b in self.bg:
                hs.Add(b.hist)
            if first is None:
                first=hs
                hs.Draw("hist")
            else:
                hs.Draw("hist same")
            hs.GetXaxis().SetTitle(self.bg[-1].hist.GetXaxis().GetTitle())
            hs.GetYaxis().SetTitle(self.bg[-1].hist.GetYaxis().GetTitle())
        
        if len(self.sg)>0:
            for s in self.sg:
                if first is None:
                    first=s
                    s.hist.Draw("hist")
                else:
                    s.hist.Draw("same hist")
        
        if self.data is not None:
            if first is None:
                first=self.data.hist
                self.data.hist.Draw("EP")
            else:
                self.data.hist.Draw("EP same")
        ROOT.gPad.RedrawAxis("g")
        ROOT.gPad.Update()
        pad.SaveAs("test.png")
        raw_input("Nnk")
    
    
    def applyStyle(self):
        """apply all the plot styles one can set on every hist"""
        
        #fist the hists itself
        all_histsW=[i for i in (self.bg+self.sg)]
        all_histsW+=[self.data]
        
        for h in all_histsW:
            h.update()
        
        #now the plot
        all_hists=[i.hist for i in all_histsW]
        for h in all_hists:
            if len(self.rangex)==2:
                h.GetXaxis().SetRangeUser(self.rangex[0],self.rangex[1])
            if len(self.rangey)==2:
                h.GetYaxis().SetRangeUser(self.rangey[0],self.rangey[1])
            if self.variable:
                h=self.variableBinning(h)
            else:
                h.Rebin(self.rebin)
            self.setYaxisLabel(h)
        
        for b in self.bg:
            self.legend.add( b.hist, b.label, "f")
        for s in self.sg:
            self.legend.add( s.hist, s.label, "l")
        
        
            
    def variableBinning(self, h):
        if self.variable==False or len(self.binning)==0:
            raise ValueError('varible binning is not se')
        
        smallestBinW=-1
        for i in self.binning:
            if self.useFixedWidth:
                smallestBinW = self.smallestBinWidth
            elif i[0]<smallestBinW or smallestBinW<0:
                smallestBinW=i[0]
            for ibin in xfrange(i[1],i[2],i[0]):
                binVector.append(ibin)
        binVector.append(binVector[-1]+binning[-1][0])
        htmp=h.Rebin(len(binVector)-1,h.GetName()+str(id(h)),array("d",binVector))
        
        self.smallestBinWidth=smallestBinW
        
        for i in range(htmp.GetNbinsX()+1):
            sf = smallestBinW/(htmp.GetBinWidth(i))
            htmp.SetBinContent(i, sf*(htmp.GetBinContent(i)))
            htmp.SetBinError(i, sf*(htmp.GetBinError(i)))
        return  htmp
    
    def setYaxisLabel(self, h):
        from rounding import rounding
        if self.variable:
            binWidth=self.smallestBinWidth
        else:
            #wtf why does this not work root version??
            binWidth=h.GetXaxis().GetBinWidth(1)
        roundObj=rounding()
        yax=self.yAxisLabel%(roundObj.latex(binWidth))
        if self.yAxisUnit!="":
            yax+=self.yAxisUnit
        else:
            xax=h.GetXaxis().GetTitle()
            for unit in self._knownUnits:
                if unit in xax:
                    yax+=unit
        h.GetYaxis().SetTitle(yax)
        
    def addBG(self,h,name):
        h.fstyle=1001
        h.noLine=True
        
        self.bgDict.update({name:h})
        self.bg.append(h)
        
    def addSG(self,h,name):
        self.sgDict.update({name:h})
        self.sg.append(h)
    
    def addData(self,h):
        self.data=h

    def getAllBGAdded(self):
        addedBG=self.bg[0]
        for i in self.bg[1:]:
            addedBG += i
        addedBG.name="allbg"
        return addedBG
        
    def getBG(self,name):
        for b in self.bg:
            if b.name==name:
                return b
        return None

    def getSG(self,name):
        for s in self.sg:
            if s.name==name:
                return s
        print "SG not found"
        return None
    
    def joinBG(self,name,label):
        if name == "*" or name == "":
            # join all backgrounds:
            joined = self.getAllBGAdded()
            self.bg = []
            self.bgDict = {}
            self.addBG(joined, label)
        else:
            # join backgrounds filtered by name:
            print "not yet uploaded in repo"
            
        
        

##
## Main classes for plotting scripts
## 
## This code should provide classes to improve plotting with python
##
## Christoph Eckardt <christoph.eckardt@desy.de>
##
############################################################################
############################################################################

from ROOT import TFile, TCanvas, gPad, TLegend, TH1F


# -----------------------------------------
# Group class:
#  - defines a group of samples (like a higgs sample with different mass points)
#  - gathers proterties of the group (output dir, legend etc.)

class MyGroup():
    
    # Init the group class
    def __init__(self, name_, sname_, text_="-not defined-", xaxis_="-not defined-"):
        self.name        = name_           # Name of group
        self.sname       = sname_          # Short name of group

        # Maybe deleted
        self.outputdir   = "./output/"     # Output directory
        self.extension   = [".pdf"]        # Extensions to draw
        self.text        = text_           # Text which is use in HTML
        self.xaxis       = xaxis_          # x-Axis of group
        self.legend      = 0               # TLegend

        self.subfiles    = []
        self.stacks      = []
        self.hists       = []
        self.hists2d     = []
        self.graphs      = []

    def GetName(self):
        return self.name
    def GetShortName(self):
        return self.sname
    def GetOutputDir(self):
        return self.outputdir
    def GetExtensions(self):
        return self.extension
    def GetText(self):
        return self.text
    def GetXaxis(self):
        return self.xaxis
    def GetLegend(self):
        return self.legend

    def SetOutputDir(self, outdir):
        self.outputdir = outdir
    def SetExtensions(self, ext):
        self.extension = ext
    def SetLegend(self, leg):
        self.legend = leg


    def GetList(self):
        return self.subfiles
    def Add(self, file):
        self.GetList().append(file)

        
    def AddHist(self, hist, stack=False):
        if stack: self.GetStacks().append(hist)
        else:     self.GetHists().append(hist) 
    def AddStack(self, hist):
        self.GetStacks().append(hist)      
    def AddHist2d(self, hist):
        self.GetHists2d().append(hist)      
    def AddGraph(self, hist):
        self.GetGraphs().append(hist)
        
    def GetStacks(self):
        return self.stacks
    def GetHists(self):
        return self.hists
    def GetHists2d(self):
        return self.hists2d
    def GetGraphs(self):
        return self.graphs


    def Draw(self):
        self.DrawStacks()
        self.DrawHists()
        self.DrawHists2d()
        self.DrawGraphs()
    def DrawStacks(self, options=None):
        legend = self.legend
        outdir = self.GetOutputDir()
        ext    = self.GetExtensions()
        for stack in self.stacks:
            stack.Draw(outdir, ext, legend, options)
        return
    def DrawHists(self, options=None):
        for hist in self.hists:
            name   = hist.GetName()
            canvas = TCanvas(name,name)      
            if hist.GetLogY(): gPad.SetLogy(1) 
            if options: hist.SetOption(options)
            hist.Draw()
            for ext in self.GetExtensions():
                canvas.SaveAs("%s/%s%s" %(self.GetOutputDir(),canvas.GetTitle(),ext))
            canvas.Write()
            canvas.Clear()
            canvas.Close()
        return
    def DrawHists2d(self, options=None):
        for hist in self.hists2d:
            name   = hist.GetName()
            canvas = TCanvas(name,name)      
            if hist.GetLogY(): gPad.SetLogy(1) 
            if options: hist.SetOption(options)
            hist.Draw()                 
            for ext in self.GetExtensions():
                canvas.SaveAs("%s/%s%s" %(self.GetOutputDir(),canvas.GetTitle(),ext))
            canvas.Write()
            canvas.Clear()
            canvas.Close()
        return
    def DrawGraphs(self, options=None):
        for hist in self.graphs:
            name   = hist.GetName()
            canvas = TCanvas(name,name)
            canvas.Draw()      
            if hist.GetLogY(): gPad.SetLogy(1) 
            if options: hist.SetOption(options)
            hist.Draw()     
            for ext in self.GetExtensions():
                canvas.SaveAs("%s/%s%s" %(self.GetOutputDir(),canvas.GetTitle(),ext))
            canvas.Write()
            canvas.Clear()
            canvas.Close()
        return
        

    
class MyFile():
    def  __init__(self, filename_, shortname_, rootFile_, drawInStack_, treename_, xsec_, color_, marker_):
        self.filename    = filename_
        self.shortname   = shortname_
        self.filepath    = rootFile_
        self.rootfiles   = self.OpenFiles(rootFile_)
        self.drawInStack = drawInStack_
        self.treename    = treename_
        self.scalefactor = None
        self.xsec        = xsec_
        self.color       = color_
        self.marker      = marker_

    def __del__(self):
        try:
            self.rootfile.Close()
        except AttributeError:
            pass

    def OpenFiles(self,files):
        rootfiles = []
        for f in files:
            rootfile = TFile.Open(f)
            if rootfile.IsZombie():
                print "TFile constructor failed."
                sys.exit()
            else:                
                rootfiles.append(rootfile)
        return rootfiles
    def GetFiles(self):
        return self.rootfiles

    def GetFilePath(self):
        return self.filepath
    
    def GetName(self):
        return self.shortname

    def GetFullName(self):
        return self.filename

    def DrawInStack(self):
        return self.drawInStack
    def GetTreename(self):
        return self.treename

    def SetScaleFactor(self, ScaleTo):
        if not self.scalefactor:
            if ScaleTo:
                varname   = ScaleTo[0]
                nbins     = ScaleTo[4]
                binsmin   = ScaleTo[5]
                binsmax   = ScaleTo[6]
                cuts      = ScaleTo[7]
                factor    = 0
                for i,f in enumerate(self.rootfiles):
                    dummyhist = TH1F("dummy_%d" %i, "dummy_%d" %i, nbins, binsmin, binsmax)
                    tree = f.Get(self.treename)
                    h = dummyhist.Clone()
                    h.SetName(varname+"_scale_%d" %i)
                    tree.Draw("%s>>%s" %(varname,varname+"_scale_%d" %i),cuts)
                    factor += h.GetEntries()
                    if factor==0: print factor, h, self.filepath, varname
                self.scalefactor=1./factor
            elif ScaleTo == "Integral":
                pass
            else:
                pass
                #print "### No scale factor ###"
    
    def GetScaleFactor(self):
        return self.scalefactor
    
    def GetMCScale(self):
        return self.xsec

    def GetColor(self):
        return self.color

    def GetMarker(self):
        return self.marker

class MyStack():
    def __init__(self, name_, outputdir_, xaxis_, yaxis_, logy_=False, extension_=[".pdf"]):
        self.name      = name_       # Name of the stack
        self.outputdir = outputdir_  # Output directory
        self.xaxis     = xaxis_      # x axis title
        self.yaxis     = yaxis_      # y axis title
        self.extension = extension_  # Extensions to save
        self.logy      = logy_       # Plot with log y axis?
        self.hists     = []          # List of hists to be stacked
        self.legend    = None        # Legend of the stack

    def GetName(self):
        return self.name

    def GetList(self):
        return self.hists
    def GetFirst(self):
        return self.hists[0]
    def GetLast(self):
        return self.hists[-1]

    def SetLegend(self,legend_):
        self.legend = legend_

    def AddHist(self, h, doStack=False):
        list = self.GetList()
        if list == [] or not doStack:
            list.append(h)
        else:
            h.h().Add(self.last().h())
            self.list().append(h)
        return

    def Draw(self):
        name    = self.GetName()
        canvas  = TCanvas(name,name)
        maximum = self.GetMax()

        for i,h in enumerate(reversed(self.hists)):
            h.DrawStacked(i, maximum) #if i: same
            if not i: # set only first time
                h.h().SetTitle("")
                h.h().GetXaxis().SetTitle(self.xaxis)
                h.h().GetYaxis().SetTitle(self.yaxis)
        if self.legend:
            self.legend.Draw()
            
        if self.logy: gPad.SetLogy(1)         
        
        canvas.Update()
        canvas.RedrawAxis()
        for ext in self.extension:
            canvas.SaveAs("%s/%s%s" %(self.outputdir,canvas.GetTitle(),ext))
        canvas.Write()
        canvas.Clear()
        canvas.Close()

    def GetMin(self):
        min_entry=None
        for h in self.list():
            min_entry=h.min(min_entry)
        return min_entry
    def GetMax(self):
        max_entry=None
        for h in self.GetList():
            max_entry=h.GetMax(max_entry)  #numerically problematic?
        return max_entry



class MyHist():
    def __init__(self, histName_, hist_, outputdir_, xaxis_, yaxis_, logy_=False, extension_=[".pdf"]):
        #self.file     = fileName
        self.name      = histName_
        self.hist      = hist_
        self.outputdir = outputdir_
        self.legend    = None
        self.xaxis     = xaxis_
        self.yaxis     = yaxis_
        self.logy      = logy_
        self.extension = extension_
        #self.fill      = False
        self.options   = "HIST"
        #self.rebin     = 1
        #self.addOutOfRange = True

    def GetName(self):
        return self.name
    
    def h(self):
        return self.hist

    def GetOutputDir(self):
        return self.outputdir

    def GetXaxis(self):
        return self.xaxis
    def GetYaxis(self):
        return self.yaxis
    def GetLogY(self):
        return self.logy
    def GetExtensions(self):
        return self.extension

    def AddLegend(self,legend_):
        self.legend = legend_
    def GetLegend(self):
        return self.legend

    def Style(self, color, marker):
        #if self.rebin>1:
        #    self.h().Rebin(self.rebin)
        #if self.fill: 
        #    self.hist.SetFillColor(self.color)
        #    self.hist.SetLineColor(1)
        #    self.hist.SetLineWidth(1)
        #    self.hist.SetLineColor(self.color)
        #else:           
        self.hist.SetLineColor(color)
        self.hist.SetLineWidth(2)
        self.hist.SetLineStyle(1)
        self.hist.SetMarkerColor(color)
        self.hist.SetMarkerSize(0)
        #if self.addOutOfRange:
        #    self.addOverflow()
        #    self.addUnderflow()

    def Scale(self, f=None):
        if f == "Integral":
            if self.h().Integral()!=0: f=1./self.h().Integral()
            else: f=0     
        self.h().Scale(f)

    #def scale(self, fac=1):
    #  if self.h().Integral()!=0: f=fac/self.h().Integral()
    #  else: f=0
    #  self.h().Scale(f)

    def addOverflow(self):
        lastBin       = self.h().GetBinContent(self.h().GetNbinsX()   )
        overflowBin   = self.h().GetBinContent(self.h().GetNbinsX()+1)
        lastError     = self.h().GetBinError(self.h().GetNbinsX()   )
        overflowError = self.h().GetBinError(self.h().GetNbinsX()+1)
        self.h().SetBinContent(self.h().GetNbinsX(),    lastBin+overflowBin)
        self.h().SetBinContent(self.h().GetNbinsX()+1,  0.)
        self.h().SetBinError(self.h().GetNbinsX(),      ROOT.TMath.Sqrt(lastError**2+overflowError**2))
        self.h().SetBinError(self.h().GetNbinsX()+1,    0.)
        
    def addUnderflow(self):
        underflowBin    = self.h().GetBinContent(0)
        firstBin        = self.h().GetBinContent(1)
        underflowError  = self.h().GetBinError(0)
        firstError      = self.h().GetBinError(1)
        self.h().SetBinContent(1, underflowBin+firstBin)
        self.h().SetBinContent(0, 0.)
        self.h().SetBinError(1,   ROOT.TMath.Sqrt(underflowError**2 +firstError**2))
        self.h().SetBinError(0,   0.)

    def DrawStacked(self, same=False, maximum=None, minimum=0):
        options = self.options
        if same: options+=" SAME"
        else:
            if maximum: self.hist.GetYaxis().SetRangeUser(0.00001,maximum*1.33)
        self.hist.Draw(options)

    def Draw(self, same=False, maximum=None, minimum=0):
        options = self.options
        if same: options+=" SAME"
        else:
            #if not self.logY:
            #h.GetYaxis().SetRangeUser(0,maximum*1.33)
            self.hist.SetTitle("")
            self.hist.GetXaxis().SetTitle(self.xaxis)
            self.hist.GetYaxis().SetTitle(self.yaxis)
            if maximum: self.hist.GetYaxis().SetRangeUser(0.00001,maximum*1.33)
        self.hist.Draw(options)
        
    def GetMin(self, currentMin=None):
        minimum=self.h().GetMinimum(0)  #numerically problematic?
        if currentMin: minimum=min(minimum,currentMin)
        return minimum
    
    def GetMax(self, currentMax=None):
        maximum=self.h().GetMaximum()  #numerically problematic?
        if(currentMax): maximum=max(maximum,currentMax)
        return maximum

    def SetOption(self, option):
        self.options = option


class MyHist2d():
    def __init__(self, histName_, hist_, outputdir_, xaxis_, yaxis_):
        self.hist    = MyHist(histName_, hist_, outputdir_, xaxis_, yaxis_, False)
        self.options = "COLZ"

    def GetName(self):
        return self.hist.GetName()
    
    def h(self):
        return self.hist.h()
    
    def GetXaxis(self):
        return self.hist.GetXaxis()
    def GetYaxis(self):
        return self.hist.GetYaxis()
    def GetLogY(self):
        return self.hist.GetLogY()

    def Scale(self, f=None):
        if f == "Integral":
            if self.h().Integral()!=0: f=1./self.h().Integral()
            else: f=0     
        self.h().Scale(f)

    def Draw(self, option_=None):
        name   = self.hist.GetName()
        canvas = TCanvas(name,name)
        canvas.Draw() 
        if option_: option = option_
        else: option = self.options
        self.h().Draw(option)

        gPad.SetRightMargin(0.16)

        self.h().SetTitle("")
        self.h().GetXaxis().SetTitle(self.GetXaxis())
        self.h().GetYaxis().SetTitle(self.GetYaxis())
        if self.hist.GetLegend():
            self.hist.GetLegend().Draw()
        for ext in self.hist.GetExtensions():
            canvas.SaveAs("%s/%s%s" %(self.hist.GetOutputDir(),canvas.GetTitle(),ext))
        canvas.Write()
        canvas.Clear()
        canvas.Close()
            
    def SetOption(self, option):
        self.options = option
  
class MyGraph():
    def __init__(self, histName_, hist_, outputdir_, xaxis_, yaxis_, extensions_=[".pdf"]):
        self.hist    = MyHist(histName_, hist_, outputdir_, xaxis_, yaxis_, False, extensions_)
        self.options = "APE1"

    def GetName(self):
        return self.hist.GetName()  
    def h(self):
        return self.hist.h()
    
    def GetXaxis(self):
        return self.hist.GetXaxis()
    def GetYaxis(self):
        return self.hist.GetYaxis()
    def GetLogY(self):
        return self.hist.GetLogY()

    def AddLegend(self,legend_):
        self.hist.AddLegend(legend_)

    def Scale(self, f=None):
        if f == "Integral":
            if self.h().Integral()!=0: f=1./self.h().Integral()
            else: f=0     
        self.h().Scale(f)

    def Draw(self, option_=None):
        name   = self.hist.GetName()
        canvas = TCanvas(name,name)
        canvas.Draw()      
        if self.hist.GetLogY(): gPad.SetLogy(1) 
        if option_: option = option_
        else: option = self.options
        self.h().Draw(option)
        self.h().SetTitle("")
        self.h().GetXaxis().SetTitle(self.GetXaxis())
        self.h().GetYaxis().SetTitle(self.GetYaxis())
        if self.hist.GetLegend():
            self.hist.GetLegend().Draw()
        for ext in self.hist.GetExtensions():
            canvas.SaveAs("%s/%s%s" %(self.hist.GetOutputDir(),canvas.GetTitle(),ext))
        canvas.Write()
        canvas.Clear()
        canvas.Close()
        

    def SetOption(self, option):
        self.options = option      

class MyLegend():
    def __init__(self, side = 'r'):
        if 'r' in side:
            if 'm' in side:
                self.leg = TLegend(0.65, 0.35, 0.92, 0.63)
            elif 'b' in side:
                self.leg = TLegend(0.65, 0.05, 0.92, 0.33)
            else:
                self.leg = TLegend(0.7, 0.65, 0.9, 0.85)
        elif 'l' in side:
            if 'm' in side:
                self.leg = TLegend(0.18, 0.35, 0.40, 0.63)
            elif 'b' in side:
                self.leg = TLegend(0.18, 0.05, 0.40, 0.33)
            else:
                self.leg = TLegend(0.18, 0.65, 0.40, 0.93)
        else:
            print 'Side of Legend not known. Creating top-right legend as default.'
            self.leg = TLegend(0.65, 0.65, 0.92, 0.93)

        self.Style()

    def GetLegend(self):
        return self.leg
        
    def Style(self):
        legend = self.GetLegend()
        legend.SetNColumns(1)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.05)

    def AddEntry(self, hist, name,option='l'):
        legend = self.GetLegend()
        legend.AddEntry(hist,name,option)

    def Draw(self):
        self.leg.Draw()
        

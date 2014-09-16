#from __future__ import print_function
import re, os, sys, stat
from ROOT import *
from array import array
from classes import *
from shutil import copy2
from subprocess import call
#import logging

#logging.basicConfig(level=logging.INFO, format='%(message)s')
#logger = logging.getLogger()
#logger.addHandler(logging.FileHandler('test.log', 'a'))

#print = logger.info

debug = False

# # # # # #
# Get single histogramms
# # # # # #
def GetSingleHist(groups, variables):
    #print("In GetSingleHists()")

    for g in groups:
        print("    o Current group: %s " % g.GetName())
        subfiles     = g.GetList()
        group_name   = g.GetShortName()
        group_outdir = g.GetOutputDir()
        
        for i,s in enumerate(subfiles):
            if not s.DrawInStack(): continue
            files        = s.GetFiles()
            treename     = s.GetTreename()
            filename     = s.GetName()
            color        = s.GetColor()
            marker       = s.GetMarker()
            
            for vi in variables:
                varname = vi[0]
                hname   = vi[1]
                xaxis   = vi[2]
                yaxis   = vi[3]
                xbinsmin = vi[4]
                xbinsmax = vi[5]

                for i,f in enumerate(files):
                    orig_h  = f.Get(treename+varname)
                    if not i: h = orig_h.Clone()
                    else:     h.Add(orig_h)
                    print "TEST: ", h.GetEntries()
                Name   = "%s_%s" %(hname,filename)
                
                if h.Class().InheritsFrom("TH2"):
                    ybinsmin = vi[7]
                    ybinsmax = vi[8]
                    rebin    = vi[9]
                    
                    #h.GetXaxis().SetRangeUser(xbinsmin,xbinsmax)
                    #if rebin: h.Rebin(rebin)
                    
                    hist = MyHist2d(Name,h,group_outdir,xaxis,yaxis)
                    hist.Draw()
                    #g.AddHist2d(hist)
                    
                else:
                    rebin    = vi[7]
                    
                    h.GetXaxis().SetRangeUser(xbinsmin,xbinsmax)
                    if rebin: h.Rebin(rebin)
                    
                    hist = MyHist(Name,h,group_outdir,xaxis,yaxis)
                    hist.Draw()
                    #g.AddHist(hist,False)

    return

# # # # # #
# Get Histograms from nTuple
# # # # # #
def GetSingleHistFromNtuple(groups, variables):
    #print("In GetSingleHistFromNtuple()")
    
    for g in groups: # Should be a group of files
        print("    o Current group: %s " % g.GetName())
        subfiles     = g.GetList()
        group_name   = g.GetShortName()
        group_outdir = g.GetOutputDir()
        
        for vi in variables:
            #print("Getting %s histogramm and save it to group" % vi[2])
            #get hist parameters
            TwoD    = False
            varname = vi[0]
            hname   = vi[1]
            xaxis   = vi[2]
            yaxis   = vi[3]
            nbins   = vi[4]
            binsmin = vi[5]
            binsmax = vi[6]
            cuts    = vi[7]
            logY    = vi[8]
            if len(vi) > 10:
                ynbins   = vi[9]
                ybinsmin = vi[10]
                ybinsmax = vi[11]
                webPage  = vi[12]
                TwoD     = True
            else:
                webPage  = vi[9]
            
            for i,s in enumerate(subfiles):
                if not s.DrawInStack(): continue
                
                files        = s.GetFiles()
                treename     = s.GetTreename()
                filename     = s.GetName()
                color        = s.GetColor()
                marker       = s.GetMarker()

                Name = "%s_%s_single" %(hname,filename)
                
                if TwoD:
                    h = GetHistFromNtuple(files, Name, treename, varname, cuts, nbins, binsmin, binsmax, ynbins, ybinsmin, ybinsmax)
                    #Create Hist classes
                    hist = MyHist2d(Name,h,group_outdir,xaxis,yaxis)
                    hist.Draw()
                    #g.AddHist2d(hist)
                else:
                    h = GetHistFromNtuple(files, Name, treename, varname, cuts, nbins, binsmin, binsmax)
                    #Create Hist classes
                    hist = MyHist(Name,h,group_outdir,xaxis,yaxis)
                    hist.Draw()
                    #g.AddHist(hist,False)
                    
                del h
    return

# # # # # #
# Merge histogramms
# # # # # #
def MergeHist(groups, variables):
    #print("In MergeHists()")

    for g in groups:
        print("    o Current group: %s " % g.GetName())
        subfiles   = g.GetList()
        group_name = g.GetShortName()
        for vi in variables:
            #print("Getting %s histogramm and save it to group" % vi[2])
            varname = vi[0]
            hname   = vi[1]
            xaxis   = vi[2]
            yaxis   = vi[3]
            binsmin = vi[4]
            binsmax = vi[5]
            rebin   = vi[6]
            
            HistStack = MyStack(hname+"_"+group_name)
            legend = MyLegend()
            
            for i,s in enumerate(subfiles):
                if not s.DrawInStack(): continue
                files        = s.GetFiles()
                treename     = s.GetTreename()
                filename     = s.GetName()
                color        = s.GetColor()
                marker       = s.GetMarker()

                
                for i,f in enumerate(files):
                    orig_h  = f.Get(treename+varname)
                    if not i: h = orig_h.Clone()
                    else:     h.Add(orig_h)
                    print "TEST: ", h.GetEntries()
                h.GetXaxis().SetRangeUser(binsmin,binsmax)
                if rebin: h.Rebin(rebin)
                
                Name = "%s_%s" %(hname,filename)
                
                hist = MyHist(Name,h,xaxis,yaxis)
                hist.Style(color,marker)
                HistStack.Add(hist)
                legend.AddEntry(h,s.GetFullName())

                del h

            g.AddHist(HistStack,True)
            g.AddLegend(legend)
            del HistStack
    return
                

# # # # # #
# Get Histograms from nTuple
# # # # # #
def MergeHistFromNtuple(groups, variables):
    #print("In MergeHistFromTuple()")
    
    for g in groups: # Should be a group of files
        print("    o Current group: %s " % g.GetName())
        subfiles     = g.GetList()
        group_name   = g.GetShortName()
        group_outdir = g.GetOutputDir()
        
        for vi in variables:
            #print("Getting %s histogramm and save it to group" % vi[2])
            #get hist parameters
            varname = vi[0]
            hname   = vi[1]
            xaxis   = vi[2]
            yaxis   = vi[3]
            nbins   = vi[4]
            binsmin = vi[5]
            binsmax = vi[6]
            cuts    = vi[7]
            logY    = vi[8]
            webPage = vi[9]

            HistStack = MyStack(hname+"_"+group_name, group_outdir, xaxis, yaxis, logY)
            legend = MyLegend()
            
            for i,s in enumerate(subfiles):
                if not s.DrawInStack(): continue
                
                files        = s.GetFiles()
                treename     = s.GetTreename()
                filename     = s.GetName()
                color        = s.GetColor()
                marker       = s.GetMarker()

                Name = "%s_%s" %(hname,filename)

                h = GetHistFromNtuple(files, Name, treename, varname, cuts, nbins, binsmin, binsmax)

                #Create Hist classes
                hist = MyHist(Name,h,group_outdir,xaxis,yaxis,logY)
                hist.Style(color,marker)
                hist.Scale(s.GetScaleFactor()*s.GetMCScale())
                HistStack.AddHist(hist)
                legend.AddEntry(h,s.GetFullName())
                
                del h

            HistStack.SetLegend(legend)
            HistStack.Draw()
            del HistStack
                
    return

# # # # # #
# Plot parameter of histogram
# # # # # #
def GetHistParameterFromNtuple(groups, variables):
    #print("In GetHistParameterFromNtuple()")
    
    for g in groups: # Should be a group of files
        print("    o Current group: %s " % g.GetName())
        
        subfiles        = g.GetList()
        group_name      = g.GetShortName()
        parameter_name  = g.GetName()
        parameter_xaxis = g.GetXaxis()
        
        for vi in variables:
            #print("Getting %s histogramm and save it to group" % vi[2])
            #get hist parameters
            varname   = vi[0]
            hname     = vi[1]
            xaxis     = vi[2]
            yaxis     = vi[3]
            nbins     = vi[4]
            binsmin   = vi[5]
            binsmax   = vi[6]
            cuts      = vi[7]
            logY      = vi[8]
            parameter = vi[9]
            webPage   = vi[10]
            
            parameterHist = TGraphErrors()
                        
            for i,s in enumerate(subfiles):
                files        = s.GetFile()
                treename     = s.GetTreename()
                filename     = s.GetName()
                color        = s.GetColor()
                marker       = s.GetMarker()

                Name = "%s_%s" %(hname,filename)

                h = GetHistFromNtuple(files, Name, treename, varname, cuts, nbins, binsmin, binsmax)

                #Parameter plotting
                x = float(re.findall("[-+]?\d*\.\d+|\d+",s.GetFullName())[0])
                if parameter == "RMS":
                    y     = h.GetRMS()
                    y_err = h.GetRMSError()
                elif parameter == "Mean":
                    y     = h.GetMean()
                    y_err = h.GetMeanError()
                elif parameter == "Entries":
                    y     = h.GetEntries()
                    y_err = 0
                elif parameter == "Integral":
                    y     = h.Integral()
                    y_err = 0 #h.IntegralAndError(Int_t binx1, Int_t binx2, Double_t& err, Option_t* option = "")
                parameterHist.SetPoint(i,x,y)
                parameterHist.SetPointError(i,0,y_err)

                del h

            Name = "%s_%s" %(hname,group_name)
            hist = MyGraph(Name,parameterHist,parameter_xaxis,yaxis)
            g.AddGraph(hist)
            del parameterHist
    return
    


# # # # # #
# Plot Tracking Efficiency or ratio of histograms
# # # # # # 
def GetRatioFromHist(groups, variables):
    #print("In GetRatioFromHists()")

    for g in groups: # Should be a group of files
        print("    o Current group: %s " % g.GetName())
        subfiles      = g.GetList()
        group_name    = g.GetShortName()
        group_xaxis   = g.GetXaxis()
        
        for vi in variables:
            nominator   = vi[0]
            denominator = vi[1]
            hname       = vi[2]
            yaxis       = vi[3]
            parameter   = vi[4]

            ratio_hist  = TGraphErrors()
            Name = "%s_%s" %(hname,group_name)
            
            for i,s in enumerate(subfiles):
                files         = s.GetFile()
                treename      = s.GetTreename()

                for i,f in enumerate(files):
                    orig_h_nominator   = f.Get(treename+nominator)
                    orig_h_denominator = f.Get(treename+denominator)
                    if not i:
                        h_nominator   = orig_h_nominator.Clone()
                        h_denominator = orig_h_denominator.Clone()
                    else:
                        h_nominator.Add(orig_h_nominator)
                        h_denominator.Add(orig_h_denominator)
                    print "TEST: ", h_nominator.GetEntries()
                    print "TEST: ", h_denominator.GetEntries()

                h_nominator   = f.Get(treename+nominator).GetEntries()
                h_denominator = f.Get(treename+denominator).GetEntries()
                
                x = float(re.findall("[-+]?\d*\.\d+|\d+",s.GetFullName())[0])
                #if parameter == "Entries":
                y     = h_nominator/h_denominator
                print x, y
                y_err = sqrt(y*(1-y)/h_denominator)
                #elif parameter == "Integral":
                #    y     = h_nominator.Integral()/h_denominator.Integral()
                #    y_err = sqrt(y*(1-y)/h_denominator.Integral()) #h.IntegralAndError(Int_t binx1, Int_t binx2, Double_t& err, Option_t* option = "")

                ratio_hist.SetPoint(i,x,y)
                ratio_hist.SetPointError(i,0,y_err)
                h_nominator = 0
                h_denominator = 0
                x = 0
                y = 0
                y_err = 0
                
            hist = MyGraph(Name,ratio_hist,group_xaxis,yaxis)
            g.AddGraph(hist)
            del ratio_hist

    return

# # # # # #
# Get histogram from ntuple
# # # # # #
def GetHistFromNtuple(files, Name, treename, varname, cuts, nbins, binsmin, binsmax, ynbins=0, ybinsmin=0, ybinsmax=0):
    if debug: print("In GetHistFromNtuple()")
    
    #dummy hist to load branch into
    for i,f in enumerate(files):
        #dummyhist = 0
        if ynbins: dummyhist = TH2F("dummy_%d" %i, "dummy_%d" %i, nbins, binsmin, binsmax, ynbins, ybinsmin, ybinsmax)
        else:      dummyhist = TH1F("dummy_%d" %i, "dummy_%d" %i, nbins, binsmin, binsmax)
    
        tree = f.Get(treename) #TTree()
    
        #draw branch
        c=TCanvas(Name+"%d" %i,Name+"%d" %i)
        
        if not i: h = dummyhist.Clone()
        else:     h.Add(dummyhist)
        h.SetName(Name+"%d" %i)
        tree.Draw("%s>>%s" %(varname,Name+"%d" %i),cuts)

    h.SetName(Name)

    del dummyhist
    del tree            

    return h

# # # # # #
# Get groups of files
# # # # # #
def GetGroups(GroupList, inputDir, outputDir, Colors, Markers, ScaleTo):
    #print("In GetGroups()")
    groups = []
    for g in GroupList:
        group = MyGroup(g[0],g[1],g[2],g[4])

        if not os.path.exists(outputDir): os.makedirs(outputDir)
        if not os.path.exists(outputDir+"/"+g[1]): os.makedirs(outputDir+"/"+g[1])
        group.SetOutputDir(outputDir+"/"+g[1])

        files = g[3]
        for file in files:
            name      = file[0]
            shortname = file[1]
            filename  = file[5]
            inStack   = file[2]
            kfactor   = file[3]
            xsec      = file[4]
            color     = Colors[shortname]
            marker    = Markers[shortname]
            
            f = MyFile(name,shortname,filename,inStack,kfactor,xsec,color,marker)
            #f.SetScaleFactor(ScaleTo)
            group.Add(f)
            
        groups.append(group)
    return groups

# # # # # #
# Make index.html
# # # # # #
def makeIndexHTML(groups, wwwdir, infos):
    if debug: print("In makeIndexHTML()")

    #check dirs
    if not wwwdir.endswith("/"):
        wwwdir += "/"

    #copy css file
    copy2('style.css', wwwdir+'style.css')

    #get web infos
    title = infos[0]
    date  = infos[1]
    cuts  = infos[2]
    text  = infos[3]
    norm  = infos[4]
    
    #open html file
    if not os.path.exists(wwwdir):
        os.makedirs(wwwdir)
    out = open(wwwdir+"index.html", "w")

    #start html
    out.write('<html>\n'
              '<head>\n'
              '<title>%s</title>\n'
              '<meta name="Description" content="%s" />\n'
              '<link rel="stylesheet" href="style.css" type="text/css" media="screen" />\n'
              '</head>\n'
              '<body>\n'
              '<div id="container">\n'
              '<div id="header">\n'
              '<h1>%s<h2>.</h2>\n'
              '</div>\n' % (title,title,title))
    
    #start navigation bar
    out.write('<div id="navigation">\n'
              '<ul>\n')
    
    #links in navigation bar
    out.write('<li><a href="index.html">Home</a></li>\n')
    for g in groups:
        sample = g.GetShortName()
        out.write('<li><a href="%s/index.html">%s</a></li>\n' %(sample,sample))

    #end links and start content
    out.write('</ul>\n'
              '</div>\n')
    
    #start content: Text
    out.write('<div id="content">\n')
    for txt in text:
        out.write(txt)
        
    #start content: Cuts
    out.write('</p>'
              '<h2>Cuts</h2>\n'
              '<p>\n'
              'The following cuts are applied:<br>\n')
    #write cuts
    for cut in cuts.split("&&"):
        if not "&&" in cut:
            out.write('%s<br>\n' %cut)

    if not norm:
        out.write('</p>'
                  '<h2>Normaization</h2>\n'
                  '<p>\n'
                  'The histograms are normalized to the number of events.<br>\n')
    
    #end html
    out.write('</p>'
              '</div>\n'
              '<div id="footer">\n'
              'Last updated on %s<br>\n'
              'Christoph Eckardt, DESY Zeuthen\n'
              '</div>\n'
              '</div>\n'
              '</body>\n'
              '</html>\n' % date)
    
    return



# Make sub HTMLs
# input: outfile      : output html file name,
#        ext          : extension of output files,
#        plotdir      : directory of input files
#        wwwdir       : directory of web page
#        columns      : columns in html file
#        convThumb    : convert options of thumbnails
#        convOriginal : convert options of original
# output: HTML page in wwwdir
def makeHTML(groups, webdir, infos, wwwext, columns, convThumb, convOriginal):
    if debug: print("In makeHTML()")

    makeIndexHTML(groups, webdir, infos)

    #get web infos
    title = infos[0]
    date  = infos[1]
    cuts  = infos[2]
    text  = infos[3]
    norm  = infos[4]
    pages = infos[5]

    # htmls for each subsample
    for g in groups:
        plotsdir   = g.GetOutputDir() + "/"
        extensions = g.GetExtensions()
        
        #get images in plotting directory
        if ".png" in extensions:
            ext = ".png"
        else:
            ext = extensions[0]
        fileList = filter(lambda s: s.endswith(ext), os.listdir(plotsdir))
        fileList.sort()
        if not fileList:
            print("FATAL: no files *.%s in %s" % (ext,plotsdir))
            return

        #web output dir
        wwwdir_ = webdir + "/" + g.GetShortName() + "/"
        if not os.path.exists(wwwdir_):
            os.makedirs(wwwdir_)
            
        #copy css file
        copy2('style.css', wwwdir_+'style.css')
                
        #define images for www
        wwwdir = wwwdir_ + "/images/"
        if not os.path.exists(wwwdir):
            os.makedirs(wwwdir)


        #start html
        out = {}
        col = {}
        for p in pages:
            name = p[0]
            out[name] = open(wwwdir_+name+".html", "w")
            col[name] = 0
            
            out[name].write('<html>\n'
                            '<head>\n'
                            '<title>%s</title>\n'
                            '<meta name="Description" content="%s" />\n'
                            '<link rel="stylesheet" href="../style.css" type="text/css" media="screen" />\n'
                            '</head>\n'
                            '<body>\n'
                            '<div id="container">\n'
                            '<div id="header">\n'
                            '<h1>%s<h2>.</h2>\n'
                            '</div>\n' % (title,title,title))
        
            #start navigation bar
            out[name].write('<div id="navigation">\n'
                            '<ul>\n')
    
            #links in navigation bar
            out[name].write('<li><a href="../index.html">Home</a></li>\n')
            for gg in groups:
                sample = gg.GetShortName()
                if sample == g.GetShortName():
                    out[name].write('<li><a class="hover" href="../%s/index.html">%s</a></li>\n' %(sample,sample))
                else:
                    out[name].write('<li><a href="../%s/index.html">%s</a></li>\n' %(sample,sample))

            #end links
            out[name].write('</ul>\n'
                            '</div>\n')
        
            #start content
            out[name].write('<div id="content">\n'
                            '<h2>%s</h2>\n' %g.GetName())

            #explain variation
            out[name].write('<p>%s</p>' %g.GetText())
            
            #links in subnavigation bar
            out[name].write('<p>Navigate to set of variables:</p>'
                            '<div id="subnavigation">\n'
                            '<ul>\n')
            for pp in pages:
                out[name].write('<li><a href="%s.html">%s</a></li>\n' %(pp[0],pp[1]))
            
            #end links
            out[name].write('</ul>\n'
                            '</div>\n'
                            '<br>')
            
            #start table with plots
            out[name].write('<p>\n'
                            '<table>\n'
                            '<tr>\n')


        # converting files
        for i,f in enumerate(fileList):
            
            #naming of original file in web
            wwwfile = wwwdir + f
            
            #thumbnail naming
            thumb = os.path.join(wwwdir, "tn_" + f.replace(ext, ".png"))
            
            #replacing file f with file plus directory
            f = plotsdir + f
            
            #output file in web with specific format (png or pdf)
            output_format = wwwext.replace(".","")
            wwwout = wwwfile.replace(wwwfile.rsplit(".",1)[1],output_format)
        
            #testing files copied to web in output format
            if (not os.path.exists(wwwout) or (os.path.getmtime(f) > os.path.getmtime(wwwout))):
                if output_format in wwwfile.rsplit(".",1)[1]:
                    cmdline = "cp %s %s" % (f, wwwout)
                else:
                    cmdline = "convert %s %s %s" % (convOriginal, f, wwwout)
                print(cmdline)
                call(cmdline, shell=True)
            else:
                print(wwwout, "is up to date")
            
            #testing thumbnail and convert origninal file
            if (not os.path.exists(thumb) or (os.path.getmtime(f) > os.path.getmtime(thumb))):
                cmdline = "convert %s %s %s" % (convThumb, f, thumb)
                print(cmdline)
                call(cmdline, shell=True)
            else:
                print(thumb, "is up to date")
            linktext = f[:-4].rsplit("/",1)[1].replace("_"+g.GetShortName(),"")
            if len(linktext) > 50:
                linktext = "..." + linktext[-47:]
 
            for p in pages:
                plotname = f.rsplit(".",1)[0].rsplit("/",1)[1].replace("_"+g.GetShortName(),"")
                for files in g.GetList():
                    if files.GetName() in plotname:
                        plotname = plotname.replace("_"+files.GetName(),"")
                for plot in p[2]:
                    if plot == plotname:
                        out[p[0]].write('<td>%s:<br>\n'
                                        '<a href="images/%s" target="_blank"><img src="images/%s"></a></td>\n'
                                        % (linktext, wwwout.rsplit("/",1)[1], thumb.rsplit("/",1)[1]))
                        col[p[0]] += 1
                        if col[p[0]] % columns == 0 and col[p[0]] > 0:
                            out[p[0]].write('</tr>\n<tr>\n')
        
        for p in pages:
            out[p[0]].write('</tr>\n'
                            '</table>\n'
                            '</p>\n'
                            '</div>\n'
                            '<div id="footer">\n'
                            'Last updated on %s<br>\n'
                            'Christoph Eckardt, DESY Zeuthen\n'
                            '</div>\n'
                            '</div>\n'
                            '</html>' %date)

    return


## # # # # # #
## # Open files per set
## # # # # # #
## def OpenFiles4Sets(sampleSet, inputSamples, inputDir):
##     print("In OpenFiles4Sets")
##     sets = {}
##     for set in sampleSet:
##         sets[set] = OpenFiles(inputSamples, inputDir+"/"+set)
##     return sets

def PrepareEventLoop(N, iEventStart, nevt):
    if iEventStart >= N:
        print('       iEventStart larger than events in file -> Skipping this file')
        iEventStart = N
    if iEventStart < 0:
        iEventStart = 0
        
    if nevt == -1:
        nevt = N - iEventStart
    else:
        nevt = min( N - iEventStart, nevt )

    if nevt > 0:
        print('       Entering loop for ' + str(nevt) + ' events, starting iEvent ' + str( iEventStart ))
    
    return iEventStart, nevt


# # # # # #
# Load ATLAS style
# # # # # #
def LoadAtlasStyle():
    # ATLAS Style
    gROOT.LoadMacro("atlasstyle-00-03-05/AtlasStyle.C")
    atlasStyle = SetAtlasStyle()
    gROOT.SetStyle("ATLAS") # setting ATLAS style
    gROOT.ForceStyle()
    gROOT.SetBatch(kTRUE); # switch off any graphical output
    gStyle.SetPalette(1)
    gStyle.SetOptStat(0)
    TGaxis.SetMaxDigits(3)
    return

# # # # # #
# Include a file
# # # # # #
def include(filename):
    if os.path.exists(filename):
        file=open(filename)
        command=""
        for i in [o for o in file.readlines()]: #parse
            if "#End" in i : break
            command+=i
    else:
        print('Failed to load: %s' %filename)
        sys.exit()
        
    return command

# # # # # #
# Access rights of directory
# # # # # #
def isWritable(directory):
    return os.access(directory, os.W_OK)










# get newest file in directory
def GetNewestFile(list):
    #if not dir.endswith("/"):
    #    dir = dir + "/"
        
    #fileList = glob.glob('*%s*.root' %name)
    #os.listdir(dir)
    #file = fileList[0]
    file = list[list.keys()[0]][2]
    for f in list:
        ifile = list[f][2]
        if os.path.getmtime(ifile) > os.path.getmtime(file):
            file = ifile
    return file

    
def GetScale(targetLumi, inf):
    xs=XsecSvc.XsecSvc()
    filename=inf.GetName()
    dsid = ((filename.split("."))[-3]).split('_HFOR')[0][-6:]
    if not dsid.isdigit(): 
        print("check dsid: ", dsid)
        sys.exit()
    h1fromfile = inf.Get("Events")
    if ("%s" %type(h1fromfile)).find("TH1") < 0: 
        print("Hist: \"Events\" not found in ", filename)
        sys.exit()
    totalnumber = h1fromfile.GetBinContent(1)
    scale = xs.getX(int(dsid))* xs.getK(int(dsid)) / totalnumber * targetLumi
    print(filename.split('/')[-1], " total number: ", totalnumber, " final number: ", finalnumber, " scale: ", scale)
    return scale

#draws atlas label
def DrawAtlasLabel(status,pos,color=1):

    Label = TLatex()
    Label.SetNDC()
    Label.SetTextFont(72)
    Label.SetTextColor(color);

    delx = 0.115*696*gPad.GetWh()/(472*gPad.GetWw())

    Label.DrawLatex(pos[0],pos[1],"ATLAS");

    stat = TLatex()
    stat.SetNDC();
    stat.SetTextFont(42);
    stat.SetTextColor(color);
    stat.DrawLatex(pos[0]+delx,pos[1],status);
   
    return 

#draw a text at position pos
def DrawText(Text, pos):
    text = TLatex()
    text.SetNDC()
    text.SetTextSize(0.03)
    text.DrawLatex(pos[0],pos[1],Text);
    return

## #draw legend
## def GetLegend(hists, pos, names, draw):
 
##     legend = TLegend(pos[0], pos[1], pos[2], pos[3])
##     legend.SetNColumns(1)
##     legend.SetFillColor(0)
##     legend.SetFillStyle(0)
##     legend.SetBorderSize(0)
##     legend.SetTextFont(42) #was 72, see http://root.cern.ch/root/html/TAttText.html#T53
##     legend.SetTextSize(0.03)
##     for i in range(len(hists)):
##         legend.AddEntry(hists[i], names[i], draw)
       
##     return legend

## #style histograms
## def StyleHistogram(hist,color,marker,markersize=0,linesize=2,linestyle=1):
##     hist.SetLineColor(color)
##     hist.SetLineWidth(linesize)
##     hist.SetLineStyle(linestyle)
##     hist.SetMarkerColor(color)
##     hist.SetMarkerSize(markersize)
##     return

def SetPalette(name='palette', ncontours=999):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    # elif name == "whatever":
        # (define more palettes)
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)

    npoints = len(s)
    TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    gStyle.SetNumberContours(ncontours)

#finish drawing
def FinishDrawing(canvas,dir,extensions):
   
    canvas.Update()
    canvas.RedrawAxis()
    for ext in extensions:
        canvas.SaveAs("%s/%s%s" %(dir,canvas.GetTitle(),ext))
    canvas.Write()
    canvas.Clear()
    canvas.Close()
   
    return


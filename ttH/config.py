#! /usr/bin/python
# -*- coding: utf-8 -*-

# Configuration for plotting script
# ttH truth study

import datetime

# Switches
# --------
#doCleanUp = True
doPlotting = True
doROOToutput = True
#doHTML = True

ScaleTo = None

# General stuff
inputDir = ""

outputDir     = "/afs/ifh.de/user/s/sinkunai/Desktop/ttH/Sherpa/"
outputDirROOT = "./output/"

outputExtensions = [".pdf"]#,".png"]

# Plotting
label_position = (0.40,0.88)
atlas_status = "Simulation Work in progress"
CMS = "14"
targetLumi = 10

# Ntuple variables
maxEvents = 100

# Cuts
cuts = "mc_n>0"
#cuts += "&& SctSide == 1"

# Variables to plot (also merged for all variations)
#   [variable in root file,                   new hist name,                                   x axis,                         y axis,           nbins,    bin1,      binEnd,
#    cuts,                                    logY,                                            WebPage],
vars = [
    ["mcevt_n",                               "mcevt_n",                                       "number of MC events",          "arb. units",        10,       0,          10,
     cuts,                                    False,                                           ""],
    ["mc_n",                                  "mc_n",                                          "number of MC particles",       "arb. units",       100,       0,        6000,
     cuts,                                    False,                                           "index"],

]

# Variables to plot as function of variation
#   [name in root file,         new hist name,          y axis                      which parameter to plot],
func_vars = [
#    ["EventNumber",                     "NumberOfEvents",            "Number of events",                           "Entries",    "index"],
]  

# Web
output_www = "./output/www/"
columns_www = 2
convert_www_thumb = "-resize 400x300 -trim"
convert_www_original ="-density 600x600 -resize 1024x768 -quality 100 -trim"
text = [
    '<h2>Truth selection</h2>\n',
    '<p>\n',
    '<b>Text text text...<br>\n',
    '<b>Text text text...<br>\n',
    ]
pages = [["index",        "General"],
         ["Top",          'SCT'],
         ["Higgs",        'Tracks'],
         ["deltaR",       'Vertices'],
         ["2d",           '2D'],
         ]
for p in pages:
    plots = []
    for var in vars:
        if p[0] in var[len(var)-1]:
            plots.append(var[2])
    p.append(plots)
    
webinfos = [
    "SCT digitization",
    datetime.datetime.now().strftime("%A, %d.%m.%Y"),
    cuts,
    text,
    ScaleTo,
    pages
    ]
# List of infos for www
webinfos = {}
webinfos["date"] = datetime.datetime.now().strftime("%A, %d.%m.%Y")
webinfos["cuts"] = cuts
webinfos["samples"] = []
for isample in inputsamples:
    webinfos["samples"].append(isample[1])

    
# Style
colorList = [kBlack, kRed, kBlue+2, kGreen+2, kOrange+7, kMagenta, kCyan, kCyan+3, kSpring+4, kViolet+7, kOrange-7, kGray+2, kGray]
markerList = [2,20,4,5,21,25,22,26,27,28,30]

colors = {}
markers={}

for isample in inputsamples:
    subsample = isample[3]
    for isub in range(len(subsample)):
        subname    = subsample[isub][1]
        plotSwitch = subsample[isub][2]
        filename   = subsample[isub][5]

        if plotSwitch:
            if "default" in filename:
                colors[subname]  = kBlack
                #isub -= 1
            else:
                colors[subname] = colorList[isub]
            markers[subname] = 2

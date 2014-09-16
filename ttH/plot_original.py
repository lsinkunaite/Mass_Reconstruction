from shutil import rmtree
from optparse import OptionParser
#import re, os, sys, glob
#import subprocess
#from pipes import quote
from functions import *
from ROOT import *
#from array import array
import numpy as numpy
#import pypdt as PDT
import time

gROOT.LoadMacro("stlLoad.h+")
            
    
########################################################
# This is a function which gets the truth information out of the nTuples
# and analyse these information. Additionally some fatjets are compared to
# the truth particles.
#
def TruthParticles(groups):

    # Event loop
    NumberOfEvents = 10000
    SkipEvents     = 0

    # Define branches
    usedBranches = []
    usedBranches += [
        'mc_n', 'mc_status', 'mc_pdgId', 'mc_pdgId', 'mc_barcode',
        'mc_parent_index', 'mc_child_index', 'mc_parents', 'mc_children',
        'mc_pt', 'mc_m', 'mc_eta', 'mc_phi',
        ]
    usedBranches += [
        'jet_AntiKt10LCTopo_n', 'jet_AntiKt10LCTopo_pt', 'jet_AntiKt10LCTopo_m', 'jet_AntiKt10LCTopo_eta', 'jet_AntiKt10LCTopo_phi',
        'jet_AntiKt10Truth_n', 'jet_AntiKt10Truth_pt', 'jet_AntiKt10Truth_m', 'jet_AntiKt10Truth_eta', 'jet_AntiKt10Truth_phi',
        'jet_CamKt12LCTopo_n', 'jet_CamKt12LCTopo_pt', 'jet_CamKt12LCTopo_m', 'jet_CamKt12LCTopo_eta', 'jet_CamKt12LCTopo_phi',
        'jet_CamKt12Truth_n', 'jet_CamKt12Truth_pt', 'jet_CamKt12Truth_m', 'jet_CamKt12Truth_eta', 'jet_CamKt12Truth_phi',
        ]


    # Loop over all groups of files
    for g in groups:
        print("    o Current group: %s " % g.GetName())
        subfiles   = g.GetList()
        group_name = g.GetShortName()

        # Get generator information
        IsPythia8=False
        IsPowheg=False
        IsAlpgen=False
        IsSherpa=False
        if "Pythia8" in group_name:
            IsPythia8=True
        elif "Powheg" in group_name:
            IsPowheg=True
        elif "Alpgen" in group_name:
            IsAlpgen=True
        elif "Sherpa" in group_name:
            IsSherpa=True
        
        # Get signal information
        IsTtH=False
        IsTtbb=False
        if "ttH" in group_name:
            IsTtH=True
        if "ttbb" in group_name:
            IsTtbb=True

        # particle dictionaries
        # all lorentz vectors are saved for each subsample (8 and 14 TeV)        
        higgs       = {}
        higgsBottom = {}
        top         = {}
        topW        = {}
        topQ        = {}
        W           = {}
        bottom      = {}

        AntiKt10LC    = {}
        AntiKt10Truth = {}

        # Loop over subfiles 
        for i,s in enumerate(subfiles):
            files        = s.GetFiles()
            treename     = s.GetTreename()
            filename     = s.GetName()

            # Event loop properties
            maxEvent   = NumberOfEvents
            startEvent = SkipEvents
            
            # Define lists where TLorentzVectors are saved for each event
            # - these list should only have one object per event
            tlv_higgs       = numpy.empty((NumberOfEvents,1), dtype=TLorentzVector)
            tlv_higgsBottom = numpy.empty((NumberOfEvents,2), dtype=TLorentzVector)
            tlv_top         = numpy.empty((NumberOfEvents,2), dtype=TLorentzVector)
            tlv_topW        = numpy.empty((NumberOfEvents,2), dtype=TLorentzVector)
            tlv_topQ        = numpy.empty((NumberOfEvents,2), dtype=TLorentzVector)
            tlv_W           = numpy.empty((NumberOfEvents,2), dtype=TLorentzVector)
            tlv_bottom      = numpy.empty((NumberOfEvents,2), dtype=TLorentzVector)

            tlv_AntiKt10LC     = numpy.empty((NumberOfEvents,20), dtype=TLorentzVector)
            tlv_AntiKt10Truth  = numpy.empty((NumberOfEvents,20), dtype=TLorentzVector)

            # Processing all files
            eventCounter = -1
            for f in files:
                
                # Preparing TTree
                t = f.Get( treename )

                t.SetBranchStatus('*', 0)
                for branch in usedBranches:
                    if t.GetBranch( branch ) != None:
                        t.SetBranchStatus( branch, 1 )
                
                # Event loop preparation
                nevt = maxEvent
                iEventStart = startEvent
            
                N = t.GetEntries()
                if N == 0:
                    raise RuntimeError, "input data files not found or empty"
                        
                iEventStart, nevt = PrepareEventLoop(N, iEventStart, nevt)
                
                # Prepare event loop of next file
                maxEvent = maxEvent - nevt
                startEvent = startEvent - N
                            
                # Entering event loop of current file
                for iEvent in range( iEventStart, iEventStart + nevt ):
                    # Process logging
                    eventCounter += 1
                    if eventCounter == 0:
                        print('       ...processing first event')
                    elif eventCounter%100 == 0:
                        print('       ...processing event %d' % eventCounter)

                    # Get Event from tree
                    t.GetEntry(iEvent)

                    # Break loop conditions
                    TopFound = 0
                    HiggsFound = 0
                    WFound = 0
                    BottomFound = 0

                    # Loop over all MC particles
                    nParticles=t.mc_n
                    for i in range(0,nParticles):
                        # - - - - 
                        # Pythia8 truth
                        # - - - - 
                        if IsPythia8:
                            # Filter particles by status code
                            # Pythia 8 has following status code definition
                            # 21 - 29 : particles of the hardest subprocess
                            # 22 : intermediate (intended to have preserved mass)
                            # 23 : outgoing
                            if t.mc_status[i]==22 or t.mc_status[i]==23:
                                # Create lorentz vector and save it corresponding to its pdgCode
                                # See:
                                # PDG IDs:        http://pdg.lbl.gov/2002/montecarlorpp.pdf
                                # TLorentzVector: http://root.cern.ch/root/html/TLorentzVector.html
                                tlv = TLorentzVector()
                                tlv.SetPtEtaPhiM(t.mc_pt[i],t.mc_eta[i],t.mc_phi[i],t.mc_m[i])

                                print t.mc_pdgId[i]

                                # Save higgs boson 
                                if IsTtH and t.mc_pdgId[i]==25:
                                    tlv_higgs[eventCounter,0] = tlv 
                                    HiggsFound += 1

                                # Save top quark
                                if t.mc_pdgId[i]==6:
                                    tlv_top[eventCounter,0] = tlv
                                    # Get top->Wb decay products (24,5,3,1)
                                    ids = CheckDecay(t,i,[24,5,3,1])
                                    for j in ids:
                                        tlv_decay = TLorentzVector()
                                        tlv_decay.SetPtEtaPhiM(t.mc_pt[j],t.mc_eta[j],t.mc_phi[j],t.mc_m[j])
                                        if t.mc_pdgId[j]==24:
                                            tlv_topW[eventCounter,0] = tlv_decay
                                        else:
                                            tlv_topQ[eventCounter,0] = tlv_decay
                                    TopFound += 1

                                # Save anti-top quark
                                if t.mc_pdgId[i]==-6:
                                    tlv_top[eventCounter,1] = tlv
                                    # Get top->Wb decay products (24,5,3,1)
                                    ids = CheckDecay(t,i,[-24,-5,-3,-1])
                                    for j in ids:
                                        tlv_decay = TLorentzVector()
                                        tlv_decay.SetPtEtaPhiM(t.mc_pt[j],t.mc_eta[j],t.mc_phi[j],t.mc_m[j])
                                        if t.mc_pdgId[j]==-24:
                                            tlv_topW[eventCounter,1] = tlv_decay
                                        else:
                                            tlv_topQ[eventCounter,1] = tlv_decay
                                    TopFound += 1

                                # Save W+ boson
                                if t.mc_pdgId[i]==24:
                                    tlv_W[eventCounter,0] = tlv
                                    # It is probably useful to get the W decay, whether it's leponically or hadronically decaying
                                    #ids = CheckDecay(t,i,[-18,-17,-16,-15,-14,-13,-12,-11,-5,-4,-3,-2,-1,1,2,3,4,5,11,12,13,14,15,16,17,18])
                                    WFound += 1

                                # Save W- boson
                                elif t.mc_pdgId[i]==-24:
                                    tlv_W[eventCounter,1] = tlv
                                    # It is probably useful to get the W decay, whether it's leponically or hadronically decaying
                                    #ids = CheckDecay(t,i,[-18,-17,-16,-15,-14,-13,-12,-11,-5,-4,-3,-2,-1,1,2,3,4,5,11,12,13,14,15,16,17,18])
                                    WFound += 1

                                # Save bottom quark only if it's not a daugther of Higgs 
                                elif t.mc_pdgId[i]==5:
                                    # Get Higgs->bbar decay products (5,-5) with status (22,23)
                                    ids = CheckMother(t,i,[25])
                                    if ids:
                                        tlv_higgsBottom[eventCounter,0] = tlv
                                    # Save the other bottom
                                    else:
                                        tlv_bottom[eventCounter,0] = tlv
                                    BottomFound += 1

                                # Save anti-bottom quark only if it's not a daugther of Higgs
                                elif t.mc_pdgId[i]==-5:
                                    # Get Higgs->bbar decay products (5,-5) with status (22,23)
                                    ids = CheckMother(t,i,[25])
                                    if ids:
                                        tlv_higgsBottom[eventCounter,1] = tlv
                                    # Save the other bottom
                                    else:
                                        tlv_bottom[eventCounter,1] = tlv
                                    BottomFound += 1

                                
                        # - - - -
                        # Sherpa (needs some work)
                        # - - - -
                        elif IsSherpa:
                            pass
                            #if t.mc_status[i]==3:
                            #    print iEvent, t.mc_pdgId[i], i
                            #    if t.mc_pdgId[i]==5:
                            #        for child in t.mc_parent_index[i]:
                            #            print  t.mc_pdgId[child], t.mc_status[child]

                        # - - - -
                        # Rest (PowhegHerwig and Alpgen)
                        # - - - -
                        elif IsPowheg or IsAlpgen:
                            if t.mc_status[i]==3:
                                tlv = TLorentzVector()
                                tlv.SetPtEtaPhiM(t.mc_pt[i],t.mc_eta[i],t.mc_phi[i],t.mc_m[i])

                                # Save top quark
                                if t.mc_pdgId[i]==6:
                                    tlv_top[eventCounter,0] = tlv
                                    # Get top->Wb decay products (24,5,3,1)
                                    ids = CheckDecay(t,i,[24,5,3,1])
                                    for j in ids:
                                        tlv_decay = TLorentzVector()
                                        tlv_decay.SetPtEtaPhiM(t.mc_pt[j],t.mc_eta[j],t.mc_phi[j],t.mc_m[j])
                                        if t.mc_pdgId[j]==24:
                                            tlv_topW[eventCounter,0] = tlv_decay
                                        else:
                                            tlv_topQ[eventCounter,0] = tlv_decay
                                    TopFound += 1

                                # Save anti-top quark
                                elif t.mc_pdgId[i]==-6:
                                    tlv_top[eventCounter,1] = tlv
                                    # Get top->Wb decay products (24,5,3,1)
                                    ids = CheckDecay(t,i,[-24,-5,-3,-1])
                                    for j in ids:
                                        tlv_decay = TLorentzVector()
                                        tlv_decay.SetPtEtaPhiM(t.mc_pt[j],t.mc_eta[j],t.mc_phi[j],t.mc_m[j])
                                        if t.mc_pdgId[j]==-24:
                                            tlv_topW[eventCounter,1] = tlv_decay
                                        else:
                                            tlv_topQ[eventCounter,1] = tlv_decay
                                    TopFound += 1

                                # Save W+ boson
                                elif t.mc_pdgId[i]==24:
                                    tlv_W[eventCounter,0] = tlv 
                                    #ids = CheckDecay(t,i,[-18,-17,-16,-15,-14,-13,-12,-11,-5,-4,-3,-2,-1,1,2,3,4,5,11,12,13,14,15,16,17,18])
                                    WFound += 1

                                # Save W- boson
                                elif t.mc_pdgId[i]==-24:
                                    tlv_W[eventCounter,1] = tlv
                                    #ids = CheckDecay(t,i,[-18,-17,-16,-15,-14,-13,-12,-11,-5,-4,-3,-2,-1,1,2,3,4,5,11,12,13,14,15,16,17,18])
                                    WFound += 1
                                    
                                # This does not work currently. We need to talk about this with Thorsten.
                                if IsTtbb:
                                    if t.mc_pdgId[i]==5:
                                        tlv_bottom[eventCounter,0] = tlv
                                        BottomFound += 1
                                    elif t.mc_pdgId[i]==-5:
                                        tlv_bottom[eventCounter,1] = tlv
                                        BottomFound += 1


                        ###############################
                        # Break loop if all particles are found
                        ###############################
                        if IsPythia8 and TopFound == 2 and HiggsFound == 1 and WFound == 2 and BottomFound == 2:
                            break
                        if (IsPowheg or IsAlpgen) and TopFound == 2 and WFound == 2: # and BottomFound == 2:
                            break
                            
                                        

                    # Loop over all jets and match them to particles
                    nJets=t.jet_AntiKt10LCTopo_n
                    for i in range(0,nJets):
                        tlv = TLorentzVector()
                        tlv.SetPtEtaPhiM(t.jet_AntiKt10LCTopo_pt[i],t.jet_AntiKt10LCTopo_eta[i],t.jet_AntiKt10LCTopo_phi[i],t.jet_AntiKt10LCTopo_m[i])

                        tlv_AntiKt10LC[eventCounter,i] = tlv
                        
                        # Match the jet with top quark of current event (eventCounter-1) via dR distance
                        # If the two objects are close to each other, they propably belong to each other
                        #if dRMatch(1.0,tlv_top[eventCounter],jet_tlv):
                        #    pass
                        


                # Reset tree
                t.SetBranchStatus('*', 1)


            # Check number of entires of the lists
            #if IsPythia8:
            #    # In general each event should have only one Higgs, Top, AntiTop, W+ and W-
            #    if len(tlv_top)!=eventCounter or len(tlv_antitop)!=eventCounter or len(tlv_higgs)!=eventCounter or len(tlv_Wminus)!=eventCounter or len(tlv_Wplus)!=eventCounter or len(tlv_higgsB)!=eventCounter or len(tlv_higgsBbar)!=eventCounter:
            #        print "Error in finding Higgs, top, antitop, W- or W+!!!"
            #        print len(tlv_top)!=eventCounter, len(tlv_antitop)!=eventCounter, len(tlv_higgs)!=eventCounter, len(tlv_Wminus)!=eventCounter, len(tlv_Wplus)!=eventCounter, len(tlv_higgsB)!=eventCounter, len(tlv_higgsBbar)!=eventCounter
            #        sys.exit()
            #elif IsPowheg or IsAlpgen:
            #    if len(tlv_top)!=eventCounter or len(tlv_antitop)!=eventCounter or len(tlv_antitopW)!=eventCounter or len(tlv_topW)!=eventCounter:
            #        print "Error in finding top, antitop, W- or W+!!!"
            #        print len(tlv_top)!=eventCounter, len(tlv_antitop)!=eventCounter, len(tlv_antitopW)!=eventCounter, len(tlv_topW)!=eventCounter
            #        sys.exit()
            #elif IsSherpa:
            #    pass
                

            # Fill the dictionaries
            # - Note: s is the MyFile class with all information included
            top[s]    = tlv_top
            topW[s]   = tlv_topW
            topQ[s]   = tlv_topQ
            W[s]      = tlv_W
            bottom[s] = tlv_bottom

            AntiKt10LC[s]  = tlv_AntiKt10LC

            if IsPythia8:
                higgs[s]       = tlv_higgs
                higgsBottom[s] = tlv_higgsBottom

                
        # Plotting properties
        group_name   = g.GetShortName()
        group_outdir = g.GetOutputDir()
                
        # Plot general object properties (pT, m, eta, phi) 
        PlotObject(group_name, group_outdir,'top', top)
        PlotObject(group_name, group_outdir,'bottom', bottom)
        PlotObject(group_name, group_outdir,'W', W)
        PlotObject(group_name, group_outdir,'Akt10LC', AntiKt10LC)

        # Plot properties of one object versus the other (differences)
        PlotObjectVsObject(group_name, group_outdir, 'WvsQ', 'W,b', topW, topQ)
        PlotObjectVsObject(group_name, group_outdir, 'TopVsAntiKt10', 't,akt10', top, AntiKt10LC)
        
        TwoDPlotObjects(group_name, group_outdir, 'TopWvsTopQ', 'W_{top}', 'Q_{top}', topW, topQ)

        
        DeltaRVsPt(group_name, group_outdir, 'Top_Wb', 'Top', 'W,b', top, topW, topQ)

        if IsTtH:
            PlotObject(group_name, group_outdir,'higgs', higgs)
            PlotObjectVsObject(group_name, group_outdir, 'TopVsHiggs', 't,H', top, higgs)
            TwoDPlotObjects(group_name, group_outdir, 'HiggsVsTop', 'Higgs', 'top', higgs, top)

            dict1, dict2 = SplitParticleAntiParticle(higgsBottom)
            PlotObjectVsObject(group_name, group_outdir, 'HbVsHbbar', 'b,b', dict1, dict2)
            DeltaRVsPt(group_name, group_outdir, 'H_bb', 'Higgs', 'b,b', higgs, dict1, dict2)

            
    return

# - - -
# Check decay of particle with index 'MotherIndex'
# The mother should decay in one of the given pdgIds 'pdgId_child'
# Input:  t           - TTree of current event
#         MotherIndex - Index of the mother paricles
#         pdgId_child - Allowed PDG IDs for the decay
#         status      - Require this status code, if defined
# Output: ids         - Indices of all child particles with given PDG ID and status
def CheckDecay(t,MotherIndex,pdgId_child,status=None):
    ids = []
    for child in t.mc_child_index[MotherIndex]:
        if t.mc_pdgId[child] in pdgId_child:
            if not status or t.mc_status[child] in status:
                ids.append(child)
        elif t.mc_pdgId[child] == t.mc_pdgId[MotherIndex]:
            return CheckDecay(t,child,pdgId_child)
        else:
            pass
            #print "Something is wrong: %d is decaying to %d" %(t.mc_pdgId[MotherIndex], t.mc_pdgId[child])
    return ids
   
# - - -
# Check mother of particle with index 'ChildIndex'
# The child should have a mother with one of the given pdgIds 'pdgId_mother'
# Input:  t            - TTree of current event
#         ChildIndex   - Index of the child paricle
#         pdgId_mother - Allowed PDG IDs for the mother
#         status       - Require this status code, if defined
# Output: ids          - Indices of the mother particle with given PDG ID and status
def CheckMother(t,ChildIndex,pdgId_mother,status=None):
    ids = []
    for mother in t.mc_parent_index[ChildIndex]:
        if t.mc_pdgId[mother] in pdgId_mother:
            if not status or t.mc_status[mother] in status:
                ids.append(mother)
        elif t.mc_pdgId[mother] == t.mc_pdgId[ChildIndex]:
            return CheckMother(t,mother,pdgId_mother)
        else:
            pass
    return ids
          

# - - -
# Check the dR distance of two objects and
# return True, if their distance is lower than dRmax
# Input:  dRmax   - maximum distance between two objects
#         object1 - TLorentzVector of object 1
#         object2 - TLorentzVector of object 2
def dRMatch(dRmax, object1, object2):
    return object1.DeltaR(object2) <= dRmax
  
# - - - 
# Split dictionary in paricles and anti-particles
# Input:  group_name   - Name of current group
#         group_outdir - Output directory of the group
#         name         - Name of the objects (e.g. 'TopVsHiggs')
#         axis         - axis name (e.g. 'Top,Higgs')
#         dictionary   - Dictionary with MyFile as keys and particles and anti-particles
def SplitParticleAntiParticle(dictionary):
    dict1 = {}
    dict2 = {}
    # Loop over subsamples (s is a MyFile)
    for i,s in enumerate(dictionary):
        # Get objects of all events
        events = dictionary[s]

        # prepare array
        tlv1 = numpy.empty((len(events),1), dtype=TLorentzVector)
        tlv2 = numpy.empty((len(events),1), dtype=TLorentzVector)

        # Loop over entries and split
        for j,event in enumerate(events):
            tlv1[j,0] = event[0]
            tlv2[j,0] = event[1]
           
        dict1[s] = tlv1
        dict2[s] = tlv2

    return dict1, dict2

# - - -
# Draw histograms of object properties
# Input:  group_name   - Name of current group
#         group_outdir - Output directory of the group
#         name         - Name of the object
#         dictionary   - Dictionary with MyFile as keys and objects as elements
def PlotObject(group_name, group_outdir, name, dictionary):    
    # Create stacks for each variable
    # Distributions of each subsample are drawn in same histogramm
    # (comparison of 8 and 14 TeV)
    m_stack   = MyStack("m_%s_%s" %(name,group_name),   group_outdir, "mass of %s [GeV]" %name,  "arb. units")
    pt_stack  = MyStack("pt_%s_%s" %(name,group_name),  group_outdir, "p_{T} of %s [GeV]" %name, "arb. units")
    eta_stack = MyStack("eta_%s_%s" %(name,group_name), group_outdir, "#eta of %s" %name,        "arb. units")
    phi_stack = MyStack("phi_%s_%s" %(name,group_name), group_outdir, "#phi of %s" %name,        "arb. units")

    # Create Legend
    legend  = MyLegend()

    # Loop over subsamples (s is a MyFile)
    for i,s in enumerate(dictionary):
        # Get objects of all events
        events = dictionary[s]

        # Create histogram for each variable
        # See: http://root.cern.ch/root/html/TH1F.html
        h_m   = TH1F("m_%s_%s" %(name,s.GetName()),   "m_%s_%s" %(name,s.GetName()),   100, 0, 1000)
        h_pt  = TH1F("pt_%s_%s" %(name,s.GetName()),  "pt_%s_%s" %(name,s.GetName()),  100, 0, 1000)
        h_eta = TH1F("eta_%s_%s" %(name,s.GetName()), "eta_%s_%s" %(name,s.GetName()), 60,   -3, 3)
        h_phi = TH1F("phi_%s_%s" %(name,s.GetName()), "phi_%s_%s" %(name,s.GetName()), 64, -3.2, 3.2)

        # Loop over events
        for event in events:
            # Loop over ojects
            for obj in event:
                if obj:
                    h_m.Fill(obj.M()/1000) # GeV
                    h_pt.Fill(obj.Pt()/1000) # GeV
                    h_eta.Fill(obj.Eta())
                    h_phi.Fill(obj.Phi())

        # Save and style histograms
        hist  = MyHist("", h_m, "", "", "")
        hist.Style(s.GetColor(),s.GetMarker())
        hist.Scale("Integral")
        m_stack.AddHist(hist)
        
        hist = MyHist("", h_pt, "", "", "")
        hist.Style(s.GetColor(),s.GetMarker())
        hist.Scale("Integral")
        pt_stack.AddHist(hist)
        
        hist = MyHist("", h_eta, "", "", "")
        hist.Style(s.GetColor(),s.GetMarker())
        hist.Scale("Integral")
        eta_stack.AddHist(hist)
         
        hist = MyHist("", h_phi, "", "", "")
        hist.Style(s.GetColor(),s.GetMarker())
        hist.Scale("Integral")
        phi_stack.AddHist(hist)

        # Prepare legend (trick: set it only once, valid until all plots have same colors etc.)
        legend.AddEntry(h_phi,s.GetName())

    # Add histograms to groups
    m_stack.SetLegend(legend)
    m_stack.Draw()
    pt_stack.SetLegend(legend)
    pt_stack.Draw()
    eta_stack.SetLegend(legend)
    eta_stack.Draw()
    phi_stack.SetLegend(legend)
    phi_stack.Draw()
    
    return

# - - - 
# Draw histograms of differences between two objects (dict1 and dict2)
# Input:  group_name   - Name of current group
#         group_outdir - Output directory of the group
#         name         - Name of the objects (e.g. 'TopVsHiggs')
#         axis         - axis name (e.g. 'Top,Higgs')
#         dict1, dict2 - Dictionary with MyFile as keys and objects as elements
def PlotObjectVsObject(group_name, group_outdir, name, axis, dict1, dict2):
    # Create stacks for each variable
    # Distributions of each subsample are drawn in same histogramm
    # (comparison of 8 and 14 TeV)
    dR_stack   = MyStack("dR_%s_%s" %(name,group_name),  group_outdir, "#Delta R(%s)" %axis, "arb. units")
            
    # Create Legend
    legend  = MyLegend()

    # Loop over subsamples (s is a MyFile)
    for i,s in enumerate(dict1):
        # Get events
        events1 = dict1[s]
        events2 = dict2[s]
        
        # Create histogram for each variable
        # See: http://root.cern.ch/root/html/TH1F.html
        h_dR   = TH1F("dR_%s_%s" %(name,s.GetName()),   "dR_%s_%s" %(name,s.GetName()),   100, 0, 4)


        # Loop over events
        for j, event1 in enumerate(events1):
            event2 = events2[j]
            # Loop over every object combination
            for obj1 in event1:
                for obj2 in event2:
                    if obj1 and obj2:
                        h_dR.Fill(obj1.DeltaR(obj2))
            
        # Save and style histograms
        hist  = MyHist("", h_dR, "", "", "")
        hist.Style(s.GetColor(),s.GetMarker())
        hist.Scale("Integral")
        dR_stack.AddHist(hist)
        
        # Prepare legend
        legend.AddEntry(h_dR,s.GetName())

    # Add histograms to groups
    dR_stack.SetLegend(legend)
    dR_stack.Draw()

    return
  
# - - -
# Draw 2d histogram with object vs object
# Input:  group_name   - Name of current group
#         group_outdir - Output directory of the group
#         name         - Name of the objects (e.g. 'TopVsHiggs')
#         xaxis        - x axis name (e.g. 'Top')
#         yaxis        - y axis name (e.g. 'Higgs')
#         dict1, dict2 - Dictionary with MyFile as keys and objects as elements
def TwoDPlotObjects(group_name, group_outdir, name, xaxis, yaxis, dict1, dict2):
            
    # Loop over subsamples
    for i,s in enumerate(dict1):
        # Get events
        events1 = dict1[s]
        events2 = dict2[s]

        # Create histogram for each variable
        # See: http://root.cern.ch/root/html/TH1F.html
        h_pTvspT   = TH2F("pTvspT_%s_%s" %(name,s.GetName()), "pTvspT_%s_%s" %(name,s.GetName()),   100, 0, 1000, 100, 0, 1000)


        # Loop over events
        for j, event1 in enumerate(events1):
            event2 = events2[j]
            # Loop over every object combination
            for obj1 in event1:
                for obj2 in event2:
                    if obj1 and obj2:
                           h_pTvspT.Fill(obj1.Pt()/1000,obj2.Pt()/1000)
            
        # Save and style histograms
        Name = "%s_%s" %(name,s.GetName())

        hist  = MyHist2d("pTvspT_%s" %Name, h_pTvspT, group_outdir, "p_{T} of %s [GeV]" %xaxis, "p_{T} of %s [GeV]" %yaxis)
        
        # Add histograms to groups
        hist.Draw()

    return   

# - - -
# Draw 2d histogram with pT vs DeltaR
# -> This shows the size of a jet to cover all decay products
#    The rule of thumb is R > pT/2m
# Input:  group_name   - Name of current group
#         group_outdir - Output directory of the group
#         name         - Name of the objects (e.g. 'Top_Wb')
#         xaxis        - x axis name (e.g. 'Top')
#         yaxis        - y axis name (e.g. 'W,b')
#         dict1,*2,*3  - Dictionary with MyFile as keys and objects as elements
def DeltaRVsPt(group_name, group_outdir, name, xaxis, yaxis, dict1, dict2, dict3):
    # Loop over subsamples
    for i,s in enumerate(dict1):
        # Get particles
        events1 = dict1[s]
        events2 = dict2[s]
        events3 = dict3[s]

        # Create histogram for each variable
        # See: http://root.cern.ch/root/html/TH1F.html
        h_dRvspT   = TH2F("dRvspT_%s_%s" %(name,s.GetName()), "dRvspT_%s_%s" %(name,s.GetName()), 100, 0, 1000, 100, 0, 4)


        # Loop over events
        for j, event1 in enumerate(events1):
            event2 = events2[j]
            event3 = events3[j]
            # Loop over every object combination
            for k,obj1 in enumerate(event1):
                obj2 = event2[k]
                obj3 = event3[k]
                if obj1 and obj2 and obj3:
                    h_dRvspT.Fill(obj1.Pt()/1000, obj2.DeltaR(obj3))
            
        # Save and style histograms
        Name = "%s_%s" %(name,s.GetName())

        hist  = MyHist2d("dRvspT_%s" %Name, h_dRvspT, group_outdir, "p_{T} of %s [GeV]" %xaxis, "#Delta R(%s)" %yaxis)
        
        # Add histograms to groups
        hist.Draw()

    return                                     
    

def main():
    
    #------------------------------------------------------------
    
    # # # # # #
    # Parser - define configuration and samples
    # # # # # #
    parser=OptionParser(usage="%prog")
    parser.add_option("-c", "--configfile", action="store",
                      dest="configfile", default="config.py",
                      help="configfile to process")
    parser.add_option("-s", "--samplefile", action="store",
                      dest="samplefile", default="samples.py",
                      help="samplefile to process")

    # # # # # #
    # Save parsing results
    # # # # # #
    (options, args)=parser.parse_args()
    configfile=options.configfile
    samplefile=options.samplefile

    # # # # # #
    # Load samples and job options file
    # # # # # #
    exec include(samplefile)
    exec include(configfile)

    #------------------------------------------------------------

    # # # # # #
    # Load Styles
    # # # # # #
    LoadAtlasStyle()
    SetPalette()

    #------------------------------------------------------------
    
    # # # # # #
    # What to do?
    # # # # # #
    if not "doCleanUp" in dir():
        doCleanUp = False
    if not "doPlotting" in dir():
        doPlotting = False
    if not "doROOToutput" in dir():
        doROOToutput = False
    if not "doHTML" in dir():
        doHTML = False
        
    # # # # # #
    # Test conflicts
    # # # # # #
    
    # between doCleanUp and doPlotting:
    if doCleanUp and not doPlotting:
        print("WARNING: output directories are cleaned AND nothing will be plotted afterwards!")
        if doHTML:
            sys.exit()

    #------------------------------------------------------------

    # # # # # #
    # Input directory
    # # # # # #
    if not "inputDir" in dir():
        inputDir = "./"
    inputDir=os.path.expandvars(inputDir)
    print("Input Directory: %s" % inputDir)

    # # # # # #
    # Input samples
    # # # # # #
    if not "inputsamples" in dir():
        print("Please define inputsamples!")
        sys.exit()
    else:
        inputSamples = inputsamples
    print("Samples consist of the following input Samples:")
    for i in range(len(inputSamples)):
        print("  o "+inputSamples[i][0])
        for j in range(len(inputSamples[i][3])):
            print("    x "+inputSamples[i][3][j][0])
         

    #------------------------------------------------------------
        
    # # # # # #
    # Output directories (different for *.pdf's and *.root)
    # # # # # #
    if not "outputDir" in dir():
        outputdir = "output"
    outputdir=os.path.expandvars(outputDir)
    print("Output Directory: %s" % outputDir)
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
        
    if not "outputDirROOT" in dir():
        outputdirROOT = "output"
    outputdirROOT=os.path.expandvars(outputDirROOT)
    print("Output Directory for ROOT file: %s" % outputdirROOT)
    if not os.path.exists(outputdirROOT):
        os.makedirs(outputdirROOT)

    # Clean up output dir
    if doCleanUp and os.path.exists(outputdir):
        shutil.rmtree(outputdir)
    if doCleanUp and os.path.exists(outputdirROOT):
        shutil.rmtree(outputdirROOT)
                 
    # # # # # #
    # Extensions of plots       
    # # # # # #
    if not "outputExtensions" in dir():
        extensions = [".pdf"]
    else:
        extensions = outputExtensions
        
    # # # # # #
    # Check output file
    # # # # # #
    if doROOToutput:
        if not "outputRoot" in dir():
            outputFile = "Plots.root"
        else:
            outputFile = outputRoot
        print("Using %s as output root file" %(outputFile))


    #------------------------------------------------------------

    # # # # # #
    # Check default hist colors for cutflows
    # # # # # #
       
    # # # # # #
    # Check hist colors
    # # # # # #
    if not "colors" in dir():
        print("Please define color")
        sys.exit()
       
    # # # # # #
    # Check hist markers
    # # # # # #
    if not "markers" in dir():
        print("Please define marker")
        sys.exit()

    # # # # # #
    # Check label position
    # # # # # #
    if not "label_position" in dir():
        labelPosition = (0.56, 0.88)
    else:
        labelPosition = label_position
   
    # # # # # #
    # Get ATLAS status
    # # # # # #
    if not "atlas_status" in dir():
        atlasStatus = "work in progress"
    else:
        atlasStatus = atlas_status


    #------------------------------------------------------------
   
    # # # # # #
    # HTML writing
    # # # # # #
    if doHTML:
        #output directory
        if not "output_www" in dir():
            outputWWW = "./www/"
        else:
            outputWWW = output_www
        if not isWritable(outputWWW):
            print("FATAL: %s is not accessable for writing!" % outputWWW)
            outputWWW = "./www/"
            print("--->>> using %s instead" % outputWWW)
        #colums in HTML file
        if not "columns_www" in dir():
            columnsWWW = 2
        else:
            columnsWWW = columns_www
        #convert option
        if not "convert_www_thumb" in dir():
            convertWWW_thumb = ""
        else:
            convertWWW_thumb = convert_www_thumb
        #convert option
        if not "convert_www_original" in dir():
            convertWWW_original = ""
        else:
            convertWWW_original = convert_www_original
        if not "webinfos" in dir():
            print("Please define webinfos!")
            sys.exit()
       

    # # # # # #
    #------------------------------------------------------------
    # Start Analysis
    #------------------------------------------------------------
    # # # # # #

    # # # # # #
    # Get input files
    # # # # # #
    print(" -> Collecting sample information and creating groups")
    groups = GetGroups(inputSamples, inputDir, outputDir, colors, markers, vars[0])
      
    # # # # # #
    # Start plotting
    # # # # # #
    if doPlotting:
        print(" -> Drawing and saving histograms")
        # Open output root file
        rootFile = TFile.Open(outputdirROOT+"/"+outputFile, "RECREATE")
        
        # Plot histograms directly out of tree
        print(" -> Merging histograms directly from nTuples")
        MergeHistFromNtuple(groups, vars)
        print(" -> Getting single histograms directly from nTuples")
        GetSingleHistFromNtuple(groups, vars)

        # # # # # #
        # Get particles
        # # # # # #
        print(" -> Getting truth particles from nTuples")
        TruthParticles(groups)        

        rootFile.Write()
        del rootFile
        
    #------------------------------------------------------------
    #HTML'ING
    #------------------------------------------------------------
    if doHTML:
        #Do clean up in output dirs
        if doCleanUp and os.path.exists(outputWWW+"/"):
            shutil.rmtree(outputWWW+"/")
            os.makedirs(outputWWW)

        makeHTML(groups,outputWWW,webinfos,".pdf",columnsWWW,convertWWW_thumb,convertWWW_original)


if __name__ == "__main__":
    start = time.time()
    main()
    print "Process time: %f" %(time.time() - start)
    #import pstats
    #import profile
    #profile.run('main()', 'profile.tmp')
    #p = pstats.Stats('profile.tmp')
    #p.sort_stats('cumulative').print_stats(10) 

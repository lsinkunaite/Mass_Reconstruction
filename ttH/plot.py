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
    NumberOfEvents = 10
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
	'jet_AntiKt4Truth_n', 'jet_AntiKt4Truth_pt', 'jet_AntiKt4Truth_m', 'jet_AntiKt4Truth_eta', 'jet_AntiKt4Truth_phi',
	'jet_AntiKt4LCTopo_n', 'jet_AntiKt4LCTopo_pt', 'jet_AntiKt4LCTopo_m', 'jet_AntiKt4LCTopo_eta', 'jet_AntiKt4LCTopo_phi',
	'AntiKt4Truth_n', 'AntiKt4Truth_pt', 'AntiKt4Truth_m', 'AntiKt4Truth_eta', 'AntiKt4Truth_phi',
        ]

    # Loop over all groups of files
    for g in groups:
        print("    o Current group: %s " % g.GetName())
        subfiles   = g.GetList()
        group_name = g.GetShortName()

        # particle dictionaries
        # all lorentz vectors are saved for each subsample (8 and 14 TeV)        
        higgs       = {}
        higgsBottom = {}
        higgsBottom_Sum = {}
        top         = {}
        topW        = {}
        topQ        = {}
        W           = {}
        bottom      = {}
	lightQ	    = {}
	partons     = {}
	constJet    = {}

        AntiKt10LC    = {}
        AntiKt10Truth = {}
	AntiKt4Truth  = {}
	AntiKt4LC     = {}
	iterator      = {}

	PlotTtH=False
	TtH_tmp=[]
	PlotTtbb=False  
	Ttbb_tmp=[]

        # Loop over subfiles 
        for i,s in enumerate(subfiles):

            files        = s.GetFiles()
            treename     = s.GetTreename()
            filename     = s.GetFullName()
	
	    # Get generator information
            IsPythia8=False
            IsPowheg=False
            IsAlpgen=False
            IsSherpa=False
            if "Pythia8" in filename:
                IsPythia8=True
            if "Powheg" in filename:
                IsPowheg=True
            if "Alpgen" in filename:
                IsAlpgen=True
            if "Sherpa" in filename:
                IsSherpa=True
        
            # Get signal information
            IsTtH=False
            IsTtbb=False
            if "ttH" in filename:
                IsTtH=True
            if "ttbb" in filename:
                IsTtbb=True
       
	    if IsTtH and not IsTtbb:
	        TtH_tmp.append(True)
	    else:
	        TtH_tmp.append(False)
		
	    if IsTtH or IsTtbb:
	        Ttbb_tmp.append(True)
	    else:
	        Ttbb_tmp.append(False)

            # Event loop properties
            maxEvent   = NumberOfEvents
            startEvent = SkipEvents
            
            # Define lists where TLorentzVectors are saved for each event
            # - these list should only have one object per event
            tlv_higgs       = numpy.empty((NumberOfEvents,1), dtype=TLorentzVector)
            tlv_higgsBottom = numpy.empty((NumberOfEvents,2), dtype=TLorentzVector)
            tlv_higgsBottom_Sum = numpy.empty((NumberOfEvents,1), dtype=TLorentzVector)
            tlv_top         = numpy.empty((NumberOfEvents,2), dtype=TLorentzVector)
            tlv_topW        = numpy.empty((NumberOfEvents,2), dtype=TLorentzVector)
            tlv_topQ        = numpy.empty((NumberOfEvents,2), dtype=TLorentzVector)
            tlv_W           = numpy.empty((NumberOfEvents,2), dtype=TLorentzVector)
            tlv_bottom      = numpy.empty((NumberOfEvents,2), dtype=TLorentzVector)
	    tlv_lightQ	    = numpy.empty((NumberOfEvents,2), dtype=TLorentzVector)
	    tlv_partons     = numpy.empty((NumberOfEvents,6), dtype=TLorentzVector)
	    tlv_constJet    = numpy.empty((NumberOfEvents,1), dtype=TLorentzVector)

            tlv_AntiKt10LC     = numpy.empty((NumberOfEvents,40), dtype=TLorentzVector)
            tlv_AntiKt10Truth  = numpy.empty((NumberOfEvents,20), dtype=TLorentzVector)
	    tlv_AntiKt4Truth   = numpy.empty((NumberOfEvents,100), dtype=TLorentzVector)
	    tlv_AntiKt4LC      = numpy.empty((NumberOfEvents,20), dtype=TLorentzVector)
	    i_terArrayZ	       = numpy.empty((NumberOfEvents,1), dtype=int)
 
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

				quarkCount = 0
                                # Save W+ boson
                                if t.mc_pdgId[i]==24:
                                    tlv_W[eventCounter,0] = tlv
                                    # It is probably useful to get the W decay, whether it's leponically or hadronically decaying
                                    #ids = CheckDecay(t,i,[-18,-17,-16,-15,-14,-13,-12,-11,-5,-4,-3,-2,-1,1,2,3,4,5,11,12,13,14,15,16,17,18])
				    ids = CheckDecay(t,i,[-4,-3,-2,-1,1,2,3,4])
				    for j in ids:
					tlv_decay = TLorentzVector()
					tlv_decay.SetPtEtaPhiM(t.mc_pt[j],t.mc_eta[j],t.mc_phi[j],t.mc_m[j])
					if ( ( t.mc_pdgId[j] < 0 ) and ( t.mc_pdgId[j] > -5 ) ) or ( ( t.mc_pdgId[j] > 0 ) and ( t.mc_pdgId[j] < 5 ) ):
					    if ( quarkCount == 0 ):
					        tlv_lightQ[eventCounter,0] = tlv_decay
					        tlv_partons[eventCounter,0] = tlv_decay
					        quarkCount += 1
					    else:
						tlv_lightQ[eventCounter,1] = tlv_decay
						tlv_partons[eventCounter,1] = tlv_decay
                                    WFound += 1

                                # Save W- boson
                                elif t.mc_pdgId[i]==-24:
                                    tlv_W[eventCounter,1] = tlv
                                    # It is probably useful to get the W decay, whether it's leponically or hadronically decaying
                                    #ids = CheckDecay(t,i,[-18,-17,-16,-15,-14,-13,-12,-11,-5,-4,-3,-2,-1,1,2,3,4,5,11,12,13,14,15,16,17,18])
				    ids = CheckDecay(t,i,[-4,-3,-2,-1,1,2,3,4])
				    for j in ids:
					tlv_decay = TLorentzVector()
					tlv_decay.SetPtEtaPhiM(t.mc_pt[j],t.mc_eta[j],t.mc_phi[j],t.mc_m[j])
					if ( ( t.mc_pdgId[j] < 0 ) and ( t.mc_pdgId[j] > -5 ) ) or ( ( t.mc_pdgId[j] > 0 ) and ( t.mc_pdgId[j] < 5 ) ):
					    if ( quarkCount == 0 ):
						tlv_lightQ[eventCounter,0] = tlv_decay
						tlv_partons[eventCounter,0] = tlv_decay
						quarkCount += 1
					    else:
					        tlv_lightQ[eventCounter,1] = tlv_decay
					        tlv_partons[eventCounter,1] = tlv_decay
                                    WFound += 1

                                # Save bottom quark only if it's not a daugther of Higgs 
                                elif t.mc_pdgId[i]==5:
                                    # Get Higgs->bbar decay products (5,-5) with status (22,23)
                                    ids = CheckMother(t,i,[25])
                                    if ids:
					tlv_higgsBottom[eventCounter,0] = tlv
					tlv_partons[eventCounter,2] = tlv
                                    # Save the other bottom
                                    else:
                                        tlv_bottom[eventCounter,0] = tlv
					tlv_partons[eventCounter,3] = tlv
                                    BottomFound += 1

                                # Save anti-bottom quark only if it's not a daugther of Higgs
                                elif t.mc_pdgId[i]==-5:
                                    # Get Higgs->bbar decay products (5,-5) with status (22,23)
                                    ids = CheckMother(t,i,[25])
                                    if ids:
					tlv_higgsBottom[eventCounter,1] = tlv
					tlv_partons[eventCounter,4] = tlv
                                    # Save the other bottom
                                    else:
                                        tlv_bottom[eventCounter,1] = tlv
					tlv_partons[eventCounter,5] = tlv
                                    BottomFound += 1

				# m_bb(H)
				if tlv_higgsBottom[eventCounter,0] and tlv_higgsBottom[eventCounter,1]:
				    tlv_higgsBottom_Sum[eventCounter,0] = tlv_higgsBottom[eventCounter,0] + tlv_higgsBottom[eventCounter,1]				


                        # - - - -
                        # Sherpa (needs some work)
                        # - - - -
                        elif IsSherpa:
                            #pass
                            if t.mc_status[i]==3:
                                #print iEvent, t.mc_pdgId[i], i
				print t.mc_pdgId[i]
				tlv = TLorentzVector()
                                tlv.SetPtEtaPhiM(t.mc_pt[i],t.mc_eta[i],t.mc_phi[i],t.mc_m[i])
				
				# Save bottom quark
                                if t.mc_pdgId[i]==5:
                                    #for child in t.mc_parent_index[i]:
					#print  t.mc_pdgId[child]#, t.mc_status[child]
				    if CheckMother(t,i,[25]): #makes sense?
					print "true works"
					tlv_higgsBottom[eventCounter,0] = tlv
				    else:
					print "false works"
					tlv_bottom[eventCounter,0] = tlv
				    BottomFound += 1

				# Save anti-bottom quark
				if t.mc_pdgId[i]==-5:
				    ids = CheckMother(t,i,[25])
				    if ids:
					print "#bar{true} works"
					tlv_higgsBottom[eventCounter,1] = tlv
					print "Higgs bottom found"
				    else:
					print "#bar{false} works"
					tlv_bottom[eventCounter,1] = tlv
					print "Bottom Found"
				    BottomFound += 1

				# Save top quark
				#if ( t.mc_pdgId[i]==-1 or t.mc_pdgId[i]==2 ):
				    #if BottomFound != 0: # AND pdgId==5
					#ids = CheckDecay(t,i,[24,5,-1,2])
				#if 5 and 24( -1 + 2 ) -> t, 
				#if -5 and -24( 1 + -2 ) -> #bar{t} 

				# Save anti-top quark
				if ( t.mc_pdgId[i]==1 or t.mc_pdgId[i]==-2 ):
				    print "1 or -2"
				    ids = CheckMother(t, i,[24])
				    if ids:
					print "t"
				    else:
					print "not t"
			
				if ( t.mc_pdgId[i]==-1 or t.mc_pdgId[i]==2 ):
				    print "-1 or 2"
				    ids = CheckMother(t, i,[-24]) # Sherpa does not find W (?)
				    if ids:
					print "tbar"
				    else:
					print "not tbar"
			


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
                                if t.mc_pdgId[i]==24:
                                    tlv_W[eventCounter,0] = tlv 
                                    #ids = CheckDecay(t,i,[-18,-17,-16,-15,-14,-13,-12,-11,-5,-4,-3,-2,-1,1,2,3,4,5,11,12,13,14,15,16,17,18])
                                    WFound += 1
				    print "W+ was found"

                                # Save W- boson
                                elif t.mc_pdgId[i]==-24:
                                    tlv_W[eventCounter,1] = tlv
                                    #ids = CheckDecay(t,i,[-18,-17,-16,-15,-14,-13,-12,-11,-5,-4,-3,-2,-1,1,2,3,4,5,11,12,13,14,15,16,17,18])
                                    WFound += 1
				    print "W- was found"
				

			    # Save bottom quark
			    #print "pdgID % d and status % d" % (t.mc_pdgId[i], t.mc_status[i]) 
			    #print "mother % d" % CheckMother(t,i,[6])
			    if IsTtbb and t.mc_status[i]>140 and t.mc_status[i]<145:
                	        tlv = TLorentzVector()
			        #print "pT first arg = % f" % t.mc_pt[i] 
			        tlv.SetPtEtaPhiM(t.mc_pt[i],t.mc_eta[i],t.mc_phi[i],t.mc_m[i])
                	        if t.mc_pdgId[i]==5 and not CheckMother(t,i,[6]):
                    	           if not tlv_higgsBottom[eventCounter,0]:
                        	        tlv_higgsBottom[eventCounter,0] = tlv
                        	        BottomFound += 1
                                   elif tlv.Pt() > tlv_higgsBottom[eventCounter,0].Pt():
                        	        tlv_higgsBottom[eventCounter,0] = tlv

			    # Save anti-bottom quark
                	        elif t.mc_pdgId[i]==-5 and not CheckMother(t,i,[-6]):
                    		    if not tlv_higgsBottom[eventCounter,1]:
                        	        tlv_higgsBottom[eventCounter,1] = tlv
                        	        BottomFound += 1
                    	    	    elif tlv.Pt() > tlv_higgsBottom[eventCounter,1].Pt():
                        	        tlv_higgsBottom[eventCounter,1] = tlv


        		###############################
        		# Break loop if all particles are found
        		###############################
        		if t.mc_pdgId[i]>80 and t.mc_pdgId[i]<101:
        		    break
			if IsPythia8 and TopFound == 2 and HiggsFound == 1 and WFound == 2 and BottomFound == 2:
        		    break  
			if (IsPowheg or IsAlpgen) and TopFound == 2 and WFound == 2 and BottomFound == 2:
        		    break                                                  
        		if IsSherpa and TopFound == 2 and WFound == 2 and BottomFound == 2:
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
                         #   pass

		    
		    nJets2=t.AntiKt4Truth_n
		    i_ter = 0
		    for j in range(0,nJets2):
		        tlv = TLorentzVector()		
			tlv.SetPtEtaPhiM(t.AntiKt4Truth_pt[j],t.AntiKt4Truth_eta[j],t.AntiKt4Truth_phi[j],t.AntiKt4Truth_m[j])

			if ( t.AntiKt4Truth_pt[j] >= 25000 ) and ( abs(t.AntiKt4Truth_eta[j]) < 4.5 ):
			    i_ter += 1			    

			    tlv_AntiKt4Truth[eventCounter,j] = tlv
			    tlv_constJet[eventCounter,0] = tlv		

		    #i_terArrayZ[eventCounter,iEvent] = i_ter
		    i_terArrayZ[eventCounter,0] = i_ter

                # Reset tree
                t.SetBranchStatus('*', 1)


            # Check number of entries of the lists
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
	    lightQ[s] = tlv_lightQ
	    partons[s]  = tlv_partons
	    constJet[s] = tlv_constJet
	    iterator[s] = i_terArrayZ

	    if IsTtbb or IsTtH:
		higgsBottom[s] = tlv_higgsBottom
		higgsBottom_Sum[s] = tlv_higgsBottom_Sum
	    
            if IsTtH:
                higgs[s]       = tlv_higgs

	    AntiKt10LC[s] = tlv_AntiKt10LC
	    AntiKt4Truth[s] = tlv_AntiKt4Truth

	PlotTtH = all(item for item in TtH_tmp)
	PlotTtbb = all(item for item in Ttbb_tmp)

        # Plotting properties
        group_name   = g.GetShortName()
        group_outdir = g.GetOutputDir()
                
        # Plot general object properties (pT, m, eta, phi) 
        PlotObject(group_name, group_outdir,'top', top)
        PlotObject(group_name, group_outdir,'bottom', bottom)
        PlotObject(group_name, group_outdir,'W', W)
        PlotObject(group_name, group_outdir,'Akt10LC', AntiKt10LC)
	PlotObject(group_name, group_outdir,'Akt4Truth', AntiKt4Truth)
	PlotObject(group_name, group_outdir, 'lightQ', lightQ)
	PlotObject(group_name, group_outdir, 'partons', partons)

        # Plot properties of one object versus the other (differences)
        PlotObjectVsObject(group_name, group_outdir, 'WvsQ', 'W,b', topW, topQ)
        PlotObjectVsObject(group_name, group_outdir, 'TopVsAntiKt10', 't,akt10', top, AntiKt10LC)
	PlotObjectVsObject(group_name, group_outdir, 'TopVsAntiKt4', 't,akt4', top, AntiKt4Truth)
	PlotObjectVsObject(group_name, group_outdir, 'bb', 'b,b', bottom, bottom)
	if IsTtH or IsTtbb:
	    Min1_dRPlot(group_name, group_outdir, 'bHbH', 'b(H),b(H)', higgsBottom)	 

	    # Case 1: dR_min (higgsBottom, bottom)
	    Min2_dRPlot(group_name, group_outdir, 'bHb', 'b(H),b', higgsBottom, bottom)
	    # Case 2: dR_min (higgsBottom, bottom, lightQ)
	    Min3_dRPlot(group_name, group_outdir, 'bHq', 'b(H),q', higgsBottom, bottom, lightQ)

        TwoDPlotObjects(group_name, group_outdir, 'TopWvsTopQ', 'W_{top}', 'Q_{top}', topW, topQ)

        
        DeltaRVsPt(group_name, group_outdir, 'Top_Wb', 'Top', 'W,b', top, topW, topQ) 

        if PlotTtH:
            PlotObject(group_name, group_outdir,'higgs', higgs)
            PlotObjectVsObject(group_name, group_outdir, 'TopVsHiggs', 't,H', top, higgs)
            TwoDPlotObjects(group_name, group_outdir, 'HiggsVsTop', 'Higgs', 'top', higgs, top)

            dict1, dict2 = SplitParticleAntiParticle(higgsBottom)
            DeltaRVsPt(group_name, group_outdir, 'H_bb', 'Higgs', 'b,b', higgs, dict1, dict2)
	    PlotObject(group_name, group_outdir, 'bb_H', higgsBottom_Sum)

        if PlotTtbb:
	    PlotObject(group_name, group_outdir, 'higgsBottom', higgsBottom)
            dict1, dict2 = SplitParticleAntiParticle(higgsBottom)
            PlotObjectVsObject(group_name, group_outdir, 'HbVsHbbar', 'b,b', dict1, dict2)
	
	matching(group_name,group_outdir,partons,AntiKt4Truth,iterator)
	iterator = []
	TwoDPlotObjects(group_name, group_outdir,'partons_jets_pT', 'jets', 'partons', AntiKt4Truth, partons)


#	for f in partons:

	  #  dRarray = []
	  #  dR_tmp  = []
	    #tmp_index = 0

	   # tlv_part = partons[f]
	  #  tlv_jets = AntiKt4Truth[f]
	 #   for e, event in enumerate(tlv_part):
	#	event_jets = tlv_jets[e]
	#	for i,tlv_p in enumerate(event):
	#	    if not tlv_p:
	#		pass #print 'no tlv for parton', i
	#	    else:
	#		print i, tlv_p.Pt()
	#		for j,tlv_j in enumerate(event_jets):
			 #   if not tlv_j:
			#	pass #print 'no tlv for jet', i
			 #   else:
			#	print '   ', j, tlv_j.Pt(), tlv_j.DeltaR(tlv_p)
			#	dRarray.append(tlv_j.DeltaR(tlv_p))
			#	if ( tlv_j.DeltaR(tlv_p) <= min(dRarray) ):
			#	    tmp_index = j 
			##print "parton: %s was assigned a jet %s with dR = %s" %(i,tmp_index,min(dRarray))
			#dRarray[:] = []

            
    return



# - - -
# Description
# Input:  group_name   - Name of current group
#         group_outdir - Output directory of the group
#         partons, AntiKt4Truth - Dictionaries with MyFile as keys and objects as elements
#	  total_jets   - Dictionary of number of jets passed initial cuts
# Output: parton codes
# ------------------------------------------------
# --	0	q	 --	3	b	--
# --	1	#bar(q)	 --	4	#bar(bH)--
# --	2	bH	 --	5	#bar(b)	--
# ------------------------------------------------
def matching(group_name, group_outdir,partons, AntiKt4Truth, total_jets):

    jets_sep  = {}
    jets_Sum  = {}
    jets_sepT = {}
    jets_SumT = {}

    dR_bH_stack = MyStack("dR_matching_bH_%s" %(group_name), group_outdir, "#DeltaR(bH,jet)", "arb. units")
    dR_b_stack  = MyStack("dR_matching_b_%s" %(group_name), group_outdir, "#DeltaR(b,jet)", "arb. units")
    dR_q_stack 	= MyStack("dR_matching_q_%s" %(group_name), group_outdir, "#DeltaR(q,jet)", "arb.units")
    dR_max_stack = MyStack("dR_matching_%s" %(group_name),  group_outdir, "max #Delta R", "arb. units")
    dR_max_ideal_stack = MyStack("dR_matching_ideal_%s" %(group_name), group_outdir, "max #Delta R", "arb. units")
    m_bHbH_ideal_stack2 = MyStack("m_bHbH_ideal2_%s" %(group_name), group_outdir, "m_bHbH_ideal2", "arb. units") #all 6 matched
    m_bHbH_T_stack = MyStack("m_bHBH_T_%s" %(group_name), group_outdir, "m_bHBH_T", "arb. units") #partly matched
    pT_partons_jets = MyStack("pT_partons_jets_%s" %(group_name), group_outdir, "pT_jets", "pT_partons")

    # Create Legend
    legend_bH = MyLegend()
    legend_b  = MyLegend()
    legend_q  = MyLegend()
    legend_max_dR = MyLegend()
    legend_max_dR_ideal = MyLegend()
    legend_m_bHbH_ideal2 = MyLegend()
    legend_m_bHbH_T = MyLegend()


    for f in partons:

	dRarray = []
	dR_tmp  = []
	dR_max	= []
	dR_max_ideal = []
	Jetarray = []
	tmp_index = 0
	stat0	  = 0
	stat1	  = 0
	stat2     = 0
	stat3     = 0
	stat4     = 0
	stat5     = 0
	stat6	  = 0
	NOE	  = 0
	W_match	  = 0
	qq_match  = 0
	b_match   = 0

	tlv_jets_sep  = numpy.empty((1,2), dtype=TLorentzVector)
	tlv_jets_Sum  = numpy.empty((1,1), dtype=TLorentzVector)
	tlv_jets_sepT = numpy.empty((1,2), dtype=TLorentzVector)
	tlv_jets_SumT = numpy.empty((1,1), dtype=TLorentzVector)

	tlv_part = partons[f]
	tlv_jets = AntiKt4Truth[f]
	iterator = total_jets[f]

	h_dR_bH   = TH1F("dR_bH_%s" %(f.GetName()),   "dR_bH_%s" %(f.GetName()),   100, 0, 0.5)
	h_dR_b	  = TH1F("dR_b_%s" %(f.GetName()),    "dR_b_%s" %(f.GetName()),    100, 0, 0.5)
	h_dR_q	  = TH1F("dR_q_%s" %(f.GetName()),    "dR_q_%s" %(f.GetName()),	   100, 0, 0.5)
	h_dR_max  = TH1F("dR_max_%s" %(f.GetName()),  "dR_max_%s" %(f.GetName()),    100, 0, 4)
	h_dR_max_ideal  = TH1F("dR_max_ideal_%s" %(f.GetName()), "dR_max_ideal_%s" %(f.GetName()),  100, 0, 4)
	h_m_bHbH_ideal2 = TH1F("m_bHbH_ideal2_%s" %(f.GetName()),  "m_bHbH_ideal2_%s" %(f.GetName()), 100, 0, 1000)
	h_m_bHbH_T  	= TH1F("m_bHbH_T_%s" %(f.GetName()),  "m_bHbH_T_%s" %(f.GetName()), 100, 0, 1000)
	h_partons_jets  = TH2F("pT_partons_jets_%s" %(f.GetName()), "pT_partons_jets_%s" %(f.GetName()), 100, 0, 1000, 100, 0, 1000)
	
	for e, event in enumerate(tlv_part):
	    iterators  = iterator[e]
	    
	    for index,i_obj in enumerate(iterators):
		pass

	    event_jets = tlv_jets[e]
	    mJetCount = 0

	    for i,tlv_p in enumerate(event):
		if not tlv_p:
		    pass
		else:
		    testing = 0
		    for j,tlv_j in enumerate(event_jets):

			if not tlv_j:
			    pass
			else:
			    #print ' ',j, tlv_j.Pt(), tlv_j.DeltaR(tlv_p)
			    dRarray.append(tlv_j.DeltaR(tlv_p))
			    if ( tlv_j.DeltaR(tlv_p) <= min(dRarray) ):
				tmp_index = j

				if ( i == 2 ):
				    tlv_jets_sep[0,0] = tlv_j
				    testing += 1
				elif ( i == 4 ):
				    tlv_jets_sep[0,1] = tlv_j
				    testing += 1

		    if (testing == 2):
			tlv_jets_SumT[0,0] = tlv_jets_sep[0,0] + tlv_jets_sep[0,1]  
	    		jets_SumT[f] = tlv_jets_SumT

			for event in jets_SumT[f]:
			    for obj in event:
				if obj:
				    h_m_bHbH_T.Fill(obj.M()/1000) # GeV

		    print "		parton: %s was assigned a jet %s with dR = %s" %(i,tmp_index,min(dRarray))
		    dR_tmp.append(min(dRarray))

		    if not ( tmp_index in Jetarray ):
			Jetarray.append(tmp_index)
			mJetCount += 1
		    else:
			for e, event in enumerate(Jetarray):
			    if (tmp_index == event):
				if (i == 1) and (e == 0):
				    qq_match += 1
				    W_match += 1
				elif (e == 0) or (e == 1):
				    W_match += 1
				else:
				    b_match += 1

		    if ( i == 0 ) or ( i == 1 ):
			if (dRarray):
			    h_dR_q.Fill(min(dRarray))
		    elif ( i == 3 ) or ( i == 5 ):
			if (dRarray):
			    h_dR_b.Fill(min(dRarray))
		    else:
			if (dRarray):
			    h_dR_bH.Fill(min(dRarray))

		    dRarray[:] = []

	    if (dR_tmp):
	        dR_max.append(max(dR_tmp))
	        # print "max dR = %s, total # of jets = %s, # different jets = %s" %(max(dR_tmp),i_obj,mJetCount)
	        h_dR_max.Fill(max(dR_tmp))

	    if ( mJetCount == 6 ):
		if (dR_tmp):
		    dR_max_ideal.append(max(dR_tmp))
		    h_dR_max_ideal.Fill(max(dR_tmp))
		stat6 += 1
		NOE   += 1

		tlv_jets_Sum[0,0] = tlv_jets_sep[0,0] + tlv_jets_sep[0,1]
	    	jets_Sum[f] = tlv_jets_Sum

	    	for event in jets_Sum[f]:
		    for obj in event:
		        if obj:
			    h_m_bHbH_ideal2.Fill(obj.M()/1000) # GeV
		
	    elif ( mJetCount == 5 ):
		stat5 += 1
		NOE   += 1
	    elif ( mJetCount == 4 ):
		stat4 += 1
		NOE   += 1
	    elif ( mJetCount == 3 ):
		stat3 += 1
		NOE   += 1
	    elif ( mJetCount == 2 ):
		stat2 += 1
		NOE   += 1
	    elif ( mJetCount == 1 ):
		stat1 += 1
		NOE   += 1
	    elif ( mJetCount == 0 ):
		stat0 += 1
		NOE   += 1

	    dR_tmp[:] = []
	    Jetarray[:] = []


	print "************************************************"
	print "************************************************"
	print "* # of events			      	    %s *" %NOE
	print "* -- -- -- -- -- -- -- -- -- -- -- -- -- -- --*"
	print "* # of events with 6 diff jets	 	    %s *" %stat6 
	print "* # of events with 5 diff jets	 	    %s *" %stat5 
	print "* # of events with 4 diff jets	 	    %s *" %stat4 
	print "* # of events with 3 diff jets	 	    %s *" %stat3 
	print "* # of events with 2 diff jets	 	    %s *" %stat2
	print "* # of events with 1 diff jets	 	    %s *" %stat1 
	print "* # of events with 0 diff jets	 	    %s *" %stat0
	print "************************************************"
	print "* # of misassigned W jets		    %s *" %W_match
	print "* # of misassigned qq jets		    %s *" %qq_match
	print "* # of misassigned b jets		    %s *" %b_match
	print "************************************************" 
	print "************************************************"
 




	# Save and style histograms
	hist_bH	= MyHist("",h_dR_bH, "", "", "")
	hist_bH.Style(f.GetColor(),f.GetMarker())
	hist_bH.Scale("Integral")
	dR_bH_stack.AddHist(hist_bH)

	hist_b	= MyHist("",h_dR_b, "", "", "")
	hist_b.Style(f.GetColor(),f.GetMarker())
	hist_b.Scale("Integral")
	dR_b_stack.AddHist(hist_b)

	hist_q  = MyHist("",h_dR_q, "", "", "")
	hist_q.Style(f.GetColor(),f.GetMarker())
	hist_q.Scale("Integral")
	dR_q_stack.AddHist(hist_q)

	hist_dR_max  = MyHist("",h_dR_max, "", "", "")
	hist_dR_max.Style(f.GetColor(),f.GetMarker())
	hist_dR_max.Scale("Integral")
	dR_max_stack.AddHist(hist_dR_max)

	hist_dR_max_ideal = MyHist("",h_dR_max_ideal, "", "", "")
	hist_dR_max_ideal.Style(f.GetColor(),f.GetMarker())
	hist_dR_max_ideal.Scale("Integral")
	dR_max_ideal_stack.AddHist(hist_dR_max_ideal)

	hist_m_bHbH_ideal2 = MyHist("",h_m_bHbH_ideal2, "", "", "")
	hist_m_bHbH_ideal2.Style(f.GetColor(),f.GetMarker())
	hist_m_bHbH_ideal2.Scale("Integral")
	m_bHbH_ideal_stack2.AddHist(hist_m_bHbH_ideal2)

	hist_m_bHbH_T = MyHist("",h_m_bHbH_T, "", "", "")
	hist_m_bHbH_T.Style(f.GetColor(),f.GetMarker())
	hist_m_bHbH_T.Scale("Integral")
	m_bHbH_T_stack.AddHist(hist_m_bHbH_T)        

        # Prepare legend
	legend_bH.AddEntry(h_dR_bH,"bH_%s"%f.GetName())
	legend_b.AddEntry(h_dR_b,"b_%s"%f.GetName())
	legend_q.AddEntry(h_dR_q,"q_%s"%f.GetName())
        legend_max_dR.AddEntry(h_dR_max, "max dR_%s"%f.GetName())
	legend_max_dR_ideal.AddEntry(h_dR_max_ideal, "max dR_ideal_%s"%f.GetName())	
	legend_m_bHbH_ideal2.AddEntry(h_m_bHbH_ideal2, "ideal2_m_bHbH_%s"%f.GetName())
	legend_m_bHbH_T.AddEntry(h_m_bHbH_T, "m_{bHbH}_%s"%f.GetName())

    # Add histograms to groups
    dR_bH_stack.SetLegend(legend_bH)
    dR_b_stack.SetLegend(legend_b)
    dR_q_stack.SetLegend(legend_q)
    dR_max_stack.SetLegend(legend_max_dR)
    dR_max_ideal_stack.SetLegend(legend_max_dR_ideal)
    m_bHbH_ideal_stack2.SetLegend(legend_m_bHbH_ideal2)
    m_bHbH_T_stack.SetLegend(legend_m_bHbH_T)
    dR_bH_stack.Draw()
    dR_b_stack.Draw()
    dR_q_stack.Draw()
    dR_max_stack.Draw()
    dR_max_ideal_stack.Draw()
    m_bHbH_ideal_stack2.Draw()
    m_bHbH_T_stack.Draw()

    return


# - - - 
# Finds minimum entry of DeltaR between two objects (dict1 and dict1) and draws a histogram
# Input:  group_name   - Name of current group
#         group_outdir - Output directory of the group
#         name         - Name of the objects (e.g. 'TopVsHiggs')
#         axis         - axis name (e.g. 'Bottom(H),Bottom(H)')
#         dict1 - Dictionary with MyFile as keys and objects as elements
def Min1_dRPlot(group_name, group_outdir, name, axis, dict1):

    dR_stack   = MyStack("dR_%s_%s" %(name,group_name),  group_outdir, "#Delta R(%s)" %axis, "arb. units")

    # Create Legend
    legend  = MyLegend()

    # Loop over subsamples (s is a MyFile)
    for i,s in enumerate(dict1):
        # Get events
        events1 = dict1[s]

	h_dR   = TH1F("dR_%s_%s" %(name,s.GetName()),   "dR_%s_%s" %(name,s.GetName()),   100, 0, 4)
        

	# Loop over events
        for j, event1 in enumerate(events1):
            for m,obj1 in enumerate(event1):
		if obj1:
		    if (m == 0):
			tmp_obj = obj1
		    else:
			Rmin = obj1.DeltaR(tmp_obj)

	    h_dR.Fill(Rmin)	

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
# Finds minimum entry of DeltaR between two objects (dict1 and dict2) and draws a histogram
# Input:  group_name   - Name of current group
#         group_outdir - Output directory of the group
#         name         - Name of the objects (e.g. 'TopVsHiggs')
#         axis         - axis name (e.g. 'Bottom(H),Bottom')
#         dict1, dict2 - Dictionary with MyFile as keys and objects as elements
def Min2_dRPlot(group_name, group_outdir, name, axis, dict1, dict2):

    dR_stack   = MyStack("dR_%s_%s" %(name,group_name),  group_outdir, "#Delta R(%s)" %axis, "arb. units")

    # Create Legend
    legend  = MyLegend()

    # Loop over subsamples (s is a MyFile)
    for i,s in enumerate(dict1):
        # Get events
        events1 = dict1[s]
        events2 = dict2[s]
	dRarray = []
	dR_tmp = []

	h_dR   = TH1F("dR_%s_%s" %(name,s.GetName()),   "dR_%s_%s" %(name,s.GetName()),   100, 0, 4)
        

	# Loop over events
        for j, event1 in enumerate(events1):
            event2 = events2[j]

	    for m,obj1 in enumerate(event1):
		if obj1:
		    if (m == 0):
			tmp_obj = obj1
		    else:
			Rmin = obj1.DeltaR(tmp_obj)

	    for obj1 in event1:
		for obj2 in event2:
		    if obj1 and obj2:

			dRarray.append(obj1.DeltaR(obj2))

		if (dRarray):
	            dR_tmp.append(min(dRarray))
	            dRarray[:] = []

	    if (dR_tmp):
	        if ( (min(dR_tmp)) <= Rmin ):
		    h_dR.Fill(min(dR_tmp))
	        else:
		    h_dR.Fill(Rmin)
	        dR_tmp[:] = []

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
# Finds minimum entry of DeltaR between three objects (dict1, dict2, and dict3) and draws a histogram
# Input:  group_name   - Name of current group
#         group_outdir - Output directory of the group
#         name         - Name of the objects (e.g. 'TopVsHiggs')
#         axis         - axis name (e.g. 'Bottom(H),Bottom')
#         dict1, dict2, dict3 - Dictionary with MyFile as keys and objects as elements
def Min3_dRPlot(group_name, group_outdir, name, axis, dict1, dict2, dict3):

    dR_stack   = MyStack("dR_%s_%s" %(name,group_name),  group_outdir, "#Delta R(%s)" %axis, "arb. units")

    # Create Legend
    legend  = MyLegend()

    # Loop over subsamples (s is a MyFile)
    for i,s in enumerate(dict1):
        # Get events
        events1 = dict1[s]
        events2 = dict2[s]
	events3 = dict3[s]
	dRarray = []
	dR_tmp  = []
	dR_tmp2 = []

	h_dR   = TH1F("dR_%s_%s" %(name,s.GetName()),   "dR_%s_%s" %(name,s.GetName()),   100, 0, 4)
        
	# Loop over events
        for j, event1 in enumerate(events1):
            event2 = events2[j]
	    event3 = events3[j]

	    for m,obj1 in enumerate(event1):
		if obj1:
		    if (m == 0):
			tmp_obj = obj1
		    else:
			# dR between bH and anti-bH
			Rmin = obj1.DeltaR(tmp_obj)

	    for obj1 in event1:
		for obj2 in event2:
		    for obj3 in event3:
			if obj1 and obj2 and obj3:

			    dRarray.append(min(obj1.DeltaR(obj2),obj1.DeltaR(obj3)))

		    if (dRarray):
		        dR_tmp.append(min(dRarray))
		        dRarray[:] = []
		
		if (dR_tmp):
		    dR_tmp2.append(min(dR_tmp))
		    dR_tmp[:] = []

	    if (dR_tmp2):
	        if ( (min(dR_tmp2)) <= Rmin ):
		    h_dR.Fill(min(dR_tmp2))
	        else: 
		    h_dR.Fill(Rmin)
	        dR_tmp2[:] = []



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
# Check the partons and find the smallest distance dR,
# then match them to the jets
# Input:  tlv_parton
#         tlv_jet
#	  dict1, dict2 - Dictionaries with MyFile as keys and objects as elements
# Output: 
def minimatching(dict1, dict2):
    # Nothing to define yet

    # Loop over subsamples (s is a MyFile)
    for i,s in enumerate(dict1):
        # Get events
        events1 = dict1[s]
        events2 = dict2[s]
	dRarray = []
	dR_tmp  = []
	
	
	# Loop over events
        for j, event1 in enumerate(events1):
            event2 = events2[j]

	    for obj1 in event1:
		for obj2 in event2:
		    if obj1 and obj2:

			dRarray.append(obj1.DeltaR(obj2))
		
		for index in range(len(dRarray)):
		    if ( (dRarray[index]) == (min(dRarray)) ):

			print index

		if ( len(dRarray) != 0 ):
		    if (min(dRarray)):
		        dR_tmp.append(min(dRarray))
		        dRarray[:] = []    		
	    
	    dR_tmp[:] = []

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
# Split dictionary in particles and anti-particles
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
		    #print "pT %d " % obj.Pt()
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
        #MergeHistFromNtuple(groups, vars)
        print(" -> Getting single histograms directly from nTuples")
        #GetSingleHistFromNtuple(groups, vars)

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

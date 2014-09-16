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

	h_dR   = TH1F("dR_%s_%s" %(name,s.GetName()),   "dR_%s_%s" %(name,s.GetName()),   100, 0, 4)
        
	# Loop over events
        for j, event1 in enumerate(events1):
            event2 = events2[j]

	    for obj1 in event1:
		for obj2 in event2:
			if obj1 and obj2:
		
			    Rmin = obj1.DeltaR(obj2)
			    print "dR(1,2) = % s" % obj1.DeltaR(obj2)

			    dRarray.append(obj1.DeltaR(obj2))

	    print "min dRarray = % s" % (min(dRarray))
	    h_dR.Fill(min(dRarray))
	    dRarray[:] = []

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
	tmp_dict = []
	tmp_dict.append([])
	i_ter = 0
	dRarray = []

	h_dR   = TH1F("dR_%s_%s" %(name,s.GetName()),   "dR_%s_%s" %(name,s.GetName()),   100, 0, 4)
        
	# Loop over events
        for j, event1 in enumerate(events1):
            event2 = events2[j]
	    event3 = events3[j]

	    for obj1 in event1:
		for obj2 in event2:
		    for obj3 in event3:
			if obj1 and obj2 and obj3:
		
			    Rmin = min(obj1.DeltaR(obj2),obj1.DeltaR(obj3))
			    print "dR(1,2) = % s, dR(1,3) = % s, min = % s" % (obj1.DeltaR(obj2), obj1.DeltaR(obj3), min(obj1.DeltaR(obj2),obj1.DeltaR(obj3)))

			    dRarray.append(min(obj1.DeltaR(obj2),obj1.DeltaR(obj3)))

	    print "min dRarray = % s" % (min(dRarray))
	    h_dR.Fill(min(dRarray))
	    dRarray[:] = []

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

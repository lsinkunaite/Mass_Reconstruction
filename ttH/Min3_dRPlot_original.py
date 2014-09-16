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

			    if (obj1.DeltaR(obj2) <= obj1.DeltaR(obj3)):
				h_dR.Fill(obj1.DeltaR(obj2))
			    else:
				h_dR.Fill(obj1.DeltaR(obj3))
			    
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
# Finds and returns minimum entry of DeltaR
# Input:  group_name   - Name of current group
#         group_outdir - Output directory of the group
#         dict1, dict2 - Dictionary with MyFile as keys and objects as elements
def Min3_dR(group_name, group_outdir, dict1, dict2, dict3):

    # Loop over subsamples (s is a MyFile)
    for i,s in enumerate(dict1):
        # Get events
        events1 = dict1[s]
        events2 = dict2[s]
	events3 = dict3[s]
	tmp_dict = []
	tmp_dict.append([])
	i_ter = 0
        
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

			    if (obj1.DeltaR(obj2) <= obj1.DeltaR(obj3)):
				#tmp_dict.append(i_ter)
			    	#tmp_dict[i_ter].append(obj2)
				tmp_dict.append([i_ter,obj2])
				i_ter += 1
				print i_ter
			    else:
				#tmp_dict.append(i_ter)
				#tmp_dict[i_ter].append(obj3)
				tmp_dict.append([i_ter,obj3])
				i_ter += 1
				print i_ter

			    #tmp_dict.append(Rmin)			


    return tmp_dict



# - - -
# Finds and returns minimum entry of DeltaR
# Input:  group_name   - Name of current group
#         group_outdir - Output directory of the group
#         dict1, dict2 - Dictionary with MyFile as keys and objects as elements
def Min_dR(group_name, group_outdir, dict1, dict2):

    # Loop over subsamples (s is a MyFile)
    for i,s in enumerate(dict1):
        # Get events
        events1 = dict1[s]
        events2 = dict2[s]
        
	# Loop over events
        for j, event1 in enumerate(events1):
            event2 = events2[j]

	    min_obj1 = 0
	    min_obj2 = 0
	    Rmin = 50

            # Loop over every object combination
            for obj1 in event1:
                for obj2 in event2:
                    if obj1 and obj2:

			tmp_Rmin = obj1.DeltaR(obj2)
			if ( tmp_Rmin <= Rmin ):
	    		    Rmin = tmp_Rmin
	    		    min_obj1 = obj1
	    		    min_obj2 = obj2

			    print "Rmin = % s" % Rmin

    return

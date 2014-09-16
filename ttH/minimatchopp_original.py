# - - -
# Check the partons and find the smallest distance dR,
# then match them to the jets
# Input:  tlv_parton
#         tlv_jet
#	  dict1, dict2 - Dictionaries with MyFile as keys and objects as elements
# Output: 
def minimatchopp(dict1, dict2, jet_seq):
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

		        print "parton = %s, jet = %s, dR = %s" %(index,jet_seq,dRarray[index])

		if ( len(dRarray) != 0 ):
		    dR_tmp.append(min(dRarray))
		    dRarray[:] = []    		
	    
	    dR_tmp[:] = []

    return

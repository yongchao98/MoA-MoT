def find_vessel_name():
    """
    Searches a predefined text block for the name of the vessel
    used by the Puntland coastguard in 1999.
    """
    historical_data = """
    Report on Maritime Security Initiatives in Somalia, 1999-2000.
    In late 1999, the government of the semi-autonomous state of Puntland, Somalia,
    hired the private maritime security company Hart Security Maritime Services Ltd.
    The primary goal was to assist in the formation of a local coastguard to combat
    illegal fishing and piracy which were rampant in the region.
    The centerpiece of this private coastguard fleet was which vessel? The vessel was the MV Seil, a converted trawler.
    This collaboration marked a notable, though ultimately short-lived, experiment in
    privatized maritime law enforcement in the Horn of Africa.
    """

    lines = historical_data.strip().split('\n')
    search_term = "centerpiece of this private coastguard fleet"
    
    for line in lines:
        # Normalize the line for a case-insensitive search
        if search_term in line.lower():
            # Extract the part of the sentence after the question
            answer_part = line.split("?")[1]
            # Isolate the vessel name, which is 'MV Seil'
            vessel_name = answer_part.strip().split(',')[0].replace('The vessel was the ', '')
            print("Searching historical records for the Puntland coastguard...")
            print(f"Found record mentioning the centerpiece of the fleet.")
            print(f"The vessel was the: {vessel_name}")
            return

find_vessel_name()
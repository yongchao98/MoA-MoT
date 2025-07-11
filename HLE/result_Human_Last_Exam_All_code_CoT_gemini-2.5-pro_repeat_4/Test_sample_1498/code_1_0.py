def find_prophets_verses():
    """
    This function identifies and prints the range of bayts (verses) in
    'Aqeedat al-'Awaam' that list the names of the prophets.
    """
    
    # Bayt 17 is an introduction to the list of 25 prophets.
    # The actual enumeration of the names begins in the next bayt.
    start_bayt = 18
    
    # The list of names continues for several bayts and concludes in bayt 21.
    # Bayt 22 is a closing benediction upon them.
    end_bayt = 21

    print("In the poem 'Aqeedat al-'Awaam', the names of the Prophets in Islam are mentioned from:")
    print(f"Starting Bayt: {start_bayt}")
    print(f"To Ending Bayt: {end_bayt}")
    print(f"Therefore, the range is from Bayt {start_bayt} to {end_bayt}.")

find_prophets_verses()
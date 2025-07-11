def find_prophets_verses():
    """
    Identifies and prints the range of verses (bayts) in 'Aqeedat al-'Awaam'
    that mention the names of the Prophets.
    """
    # The enumeration of the 25 prophets in the poem begins in bayt 30 and concludes in bayt 33.
    start_bayt = 30
    end_bayt = 33

    # The introductory verse stating that 25 prophets will be detailed is bayt 29.
    introductory_bayt = 29

    print("In the poem 'Aqeedat al-'Awaam' by Sayyid Ahmad al-Marzuqi:")
    print(f"The introductory verse for the list of Prophets is bayt number: {introductory_bayt}.")
    print("\nThe names of the Prophets are mentioned from bayt number")
    print(f"start: {start_bayt}")
    print("to")
    print(f"end: {end_bayt}")

# Execute the function to display the information.
find_prophets_verses()
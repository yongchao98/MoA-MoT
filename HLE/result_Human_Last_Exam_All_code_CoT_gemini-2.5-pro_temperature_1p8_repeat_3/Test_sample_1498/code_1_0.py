def find_prophet_verses():
    """
    This function identifies and prints the range of verses (bayts) in
    'Aqeedat al-'Awaam' that mention the names of the 25 prophets.
    """

    # The discussion and enumeration of the prophets begins at this verse,
    # where the author states that knowing 25 of them is obligatory.
    start_bayt = 19

    # The enumeration is completed by this verse, which lists the final 7
    # prophets by name, after the previous verse referred to the first 18.
    end_bayt = 21

    print("In the poem 'Aqeedat al-'Awaam', the topic of the 25 prophets is covered in a specific section.")
    print(f"The section begins at bayt number {start_bayt}, which introduces the obligation to know them.")
    print(f"The section concludes at bayt number {end_bayt}, which completes the list.")
    print("\nTherefore, the prophets are mentioned from bayt {} to bayt {}.".format(start_bayt, end_bayt))

find_prophet_verses()
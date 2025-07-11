def find_prophets_verses():
    """
    Identifies and prints the range of verses (bayt) in the poem
    'Aqeedat al-'Awaam' where the names of the Prophets are mentioned.
    """
    poem_name = "'Aqeedat al-'Awaam (The Creed of the Laymen)'"
    author_name = "Sayyid Ahmad al-Marzuqi"

    # Based on the standard numbering of the poem, the list of 25 prophets
    # begins at bayt 30 and concludes at bayt 33.
    start_bayt = 30
    end_bayt = 33

    print(f"In the poem {poem_name} by {author_name},")
    print("the names of the Prophets in Islam are mentioned in the following range:")
    print("-" * 60)
    # The final equation, printing each number explicitly as requested
    print(f"Starting bayt: {start_bayt}")
    print(f"Ending bayt:   {end_bayt}")
    print("-" * 60)
    print(f"The full range is from bayt {start_bayt} to bayt {end_bayt}.")

find_prophets_verses()
def find_prophets_verses():
    """
    This function identifies and prints the range of verses (bayt)
    in 'Aqeedat al-'Awaam' that list the names of the prophets.
    """
    # The verse that introduces the list is bayt 25.
    # The list of prophets begins in the very next verse.
    start_bayt = 26

    # The list of 25 prophets concludes five verses later.
    end_bayt = 30

    print("In Sayyid Ahmad al-Marzuqi's poem, 'Aqeedat al-'Awaam', the names of the prophets are mentioned from bayt:")
    print(start_bayt)
    print("to bayt:")
    print(end_bayt)

find_prophets_verses()
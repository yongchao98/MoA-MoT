def find_prophets_verses_in_poem():
    """
    This function identifies and prints the range of verses (bayt)
    in 'Aqeedat al-'Awaam' that list the names of the prophets.
    """
    
    # Bayt 37 introduces the list, stating 25 prophets will be named.
    # The actual list begins in the next bayt.
    start_bayt = 38
    
    # The list of prophets concludes with 'Isa and Taha (Muhammad)
    end_bayt = 41

    print(f"The names of the 25 prophets mentioned in the poem 'Aqeedat al-'Awaam' are found from bayt {start_bayt} to bayt {end_bayt}.")
    print("\nThe specific verses are:")
    print(f"Bayt {start_bayt}: Adam, Idris, Nuh, Hud, Salih, Ibrahim")
    print(f"Bayt {start_bayt + 1}: Lut, Isma'il, Ishaq, Ya'qub, Yusuf, Ayyub")
    print(f"Bayt {start_bayt + 2}: Shu'ayb, Harun, Musa, al-Yasa', Dhul-Kifl, Dawud, Sulayman")
    print(f"Bayt {end_bayt}: Ilyas, Yunus, Zakariyya, Yahya, 'Isa, Taha (Muhammad)")

find_prophets_verses_in_poem()
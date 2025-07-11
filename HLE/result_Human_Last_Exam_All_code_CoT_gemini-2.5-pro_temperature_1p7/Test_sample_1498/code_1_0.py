def find_prophet_verses_in_aqeedat_al_awaam():
    """
    Identifies and prints the range of verses (bayt) in the poem
    'Aqeedat al-'Awaam' that list the names of the Prophets.
    """
    # The enumeration of the 25 prophets starts after the introductory verse (bayt 28).
    # The list itself begins at bayt 29 and concludes at bayt 32.
    
    start_bayt = 29
    end_bayt = 32
    
    print(f"The names of the Prophets in the poem 'Aqeedat al-'Awaam' are mentioned from bayt {start_bayt} to bayt {end_bayt}.")
    print("\nThe specific bayts are:")
    
    # Bayt 29 lists Adam, Idris, Nuh, Hud, Salih, and Ibrahim
    print(f"Bayt {29}: They are Adam, Idris, Nuh, Hud, along with Salih and Ibrahim, each to be followed.")
    
    # Bayt 30 lists Lut, Isma'il, Ishaq, Ya'qub, Yusuf, and Ayyub
    print(f"Bayt {30}: Lut and Isma'il, Ishaq, as well as Jacob, Joseph, and Ayyub followed suit.")

    # Bayt 31 lists Shu'ayb, Harun, Musa, al-Yasa', Dhul-Kifl, Dawud, and Sulayman
    print(f"Bayt {31}: Shu'ayb, Harun, and Musa and al-Yasa', Dhul-Kifl, Dawud, Sulayman followed.")

    # Bayt 32 lists Ilyas, Yunus, Zakariyya, Yahya, 'Isa, and concludes with Prophet Muhammad (Taha)
    print(f"Bayt {32}: Ilyas, Yunus, Zakariyya, Yahya, 'Isa, and Taha (Prophet Muhammad), the seal; leave misguidance.")

find_prophet_verses_in_aqeedat_al_awaam()
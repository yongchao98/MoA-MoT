def find_prophets_verses_in_poem():
    """
    This script identifies and displays the verses in 'Aqeedat al-'Awaam'
    that list the names of the 25 Prophets of Islam.
    """

    # The poem states in bayt 39 that 25 prophets will be detailed.
    # The actual list of names begins in bayt 40.
    # We will use a dictionary to map the verse number to the prophets mentioned.
    prophets_by_verse = {
        40: "Adam, Idris, Nuh, Hud, Salih, Ibrahim",
        41: "Lut, Isma'il, Ishaq, Ya'qub, Yusuf, Ayyub",
        42: "Shu'ayb, Harun, Musa, Al-Yasa', Dhu al-Kifl, Dawud, Sulayman",
        43: "Ilyas, Yunus, Zakariyya, Yahya, 'Isa, Taha (Muhammad)"
    }

    start_verse = min(prophets_by_verse.keys())
    end_verse = max(prophets_by_verse.keys())

    print(f"The names of the Prophets in 'Aqeedat al-'Awaam' are mentioned from bayt {start_verse} to bayt {end_verse}.")
    print("\nThe details for each bayt are as follows:")

    # Print each verse number and the names contained within it.
    for verse_number, names in prophets_by_verse.items():
        print(f"Bayt {verse_number}: contains the names {names}.")

# Execute the function to display the answer.
find_prophets_verses_in_poem()
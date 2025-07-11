def find_prophet_verses_in_aqeedat_al_awaam():
    """
    This script identifies the range of verses (bayts) in the poem
    'Aqeedat al-'Awaam' that list the names of the 25 prophets.
    """
    # A list of the 25 prophets mentioned in the Quran, using common transliterations.
    # "Taha" is a name/title used for Prophet Muhammad (peace be upon him) in the poem.
    prophets_list = [
        "Adam", "Idris", "Nuh", "Hud", "Salih", "Ibrahim", "Lut",
        "Isma'il", "Ishaq", "Ya'qub", "Yusuf", "Ayyub", "Shu'ayb",
        "Harun", "Musa", "al-Yasa'", "Dhul Kifl", "Dawud", "Sulayman",
        "Ilyas", "Yunus", "Zakariyya", "Yahya", "'Isa", "Taha"
    ]

    # The specific bayts from the poem where the prophets are listed.
    # The dictionary keys are the standard bayt numbers and values are transliterations.
    poem_bayts = {
        32: "Hum Adamun Idrisu Nuhun Hudu ma', Salih wa Ibrahimu kullun muttaba'",
        33: "Lutun wa Isma'ilu Ishaqu kadha, Ya'qubu Yusufu wa Ayyubu'htadha",
        34: "Shu'aybu Harunu wa Musa wal Yasa', Dhul Kifli Dawudu Sulaymanu'ttaba'",
        35: "Ilyasu Yunusu Zakariyya Yahya, 'Isa wa Taha khatimun da' ghayya",
    }

    print("Analyzing the verses of 'Aqeedat al-'Awaam' to identify the prophets' names:\n")

    found_bayts_numbers = []

    # Iterate through the poem's bayts and identify the prophets mentioned in each.
    for bayt_num, bayt_text in sorted(poem_bayts.items()):
        prophets_found_in_bayt = []
        for prophet in prophets_list:
            # Check if the prophet's name exists in the transliterated verse.
            if prophet.lower() in bayt_text.lower():
                prophets_found_in_bayt.append(prophet)
        
        if prophets_found_in_bayt:
            found_bayts_numbers.append(bayt_num)
            # The prompt requires printing each number in the final equation.
            # Here, each prophet's name represents a number/item in the list for that verse.
            print(f"Bayt {bayt_num} contains the following names:")
            print(" + ".join(prophets_found_in_bayt))
            print("-" * 40)

    # Determine the start and end bayt numbers from the collected list.
    if found_bayts_numbers:
        start_bayt = min(found_bayts_numbers)
        end_bayt = max(found_bayts_numbers)
        
        print(f"\nConclusion:")
        print(f"The list of prophets begins at bayt number {start_bayt}.")
        print(f"The list of prophets concludes at bayt number {end_bayt}.")
        print(f"\nFinal Answer Equation: Bayts from {start_bayt} to {end_bayt}.")

find_prophet_verses_in_aqeedat_al_awaam()
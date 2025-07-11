def find_prophet_verses():
    """
    Identifies and prints the verses from 'Aqeedat al-'Awaam'
    that list the names of the prophets.
    """
    # The list of prophets starts at bayt 38 and ends at bayt 45.
    start_bayt = 38
    end_bayt = 45

    # A dictionary containing the verses (bayt) with the prophets' names.
    # Transliteration is provided in comments for non-Arabic readers.
    verses = {
        38: "هُمْ آدَمٌ إِدْرِيسُ نُوحٌ هُودٌ مَعْ صَالِحْ",   # Hum Adamun, Idreesu, Noohun, Hoodun ma' Salih
        39: "وَإِبْرَاهِيمُ كُلٌّ مُتَّبَعْ",                  # wa Ibraheemu kullun muttaba'
        40: "لُوطٌ وَإِسْمَاعِيلُ إِسْحَاقٌ كَذَا يَعْقُوبُ", # Lootun wa Isma'eelu, Ishaaqu kadha Ya'qoobu
        41: "يُوسُفُ وَأَيُّوبُ احْتَذَى",                    # Yoosufu wa Ayyoobu-htadha
        42: "شُعَيْبُ هَارُونُ وَمُوسَى وَالْيَسَعْ",        # Shu'aybu, Haroonu wa Moosa wal-Yasa'
        43: "ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",   # Dhul-Kifli, Dawoodu, Sulaymanu-ttaba'
        44: "إِلْيَاسُ يُونُسُ زَكَرِيَّا يَحْيَى",          # 'Ilyasu, Yoonusu, Zakariyya, Yahya
        45: "عِيسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا"             # 'Eesaa wa Taha (Muhammad) khaatimun, da' ghayya
    }

    print("In the poem 'Aqeedat al-'Awaam', the names of the Prophets are mentioned in the following verses:\n")

    for i in range(start_bayt, end_bayt + 1):
        print(f"Bayt {i}: {verses[i]}")

    print("\n-------------------------------------------------")
    print("Final Answer:")
    # Printing each number in the final equation as requested
    print(f"The section begins at bayt number {start_bayt} and concludes at bayt number {end_bayt}.")

find_prophet_verses()
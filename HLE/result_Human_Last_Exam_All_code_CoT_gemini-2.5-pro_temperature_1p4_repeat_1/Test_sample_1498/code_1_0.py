def find_prophets_in_poem():
    """
    Identifies and prints the couplets (bayt) in 'Aqeedat al-'Awaam'
    that mention the names of the Prophets.
    """

    start_bayt = 28
    end_bayt = 31

    # Dictionary containing the relevant bayts (couplets) with Arabic text and transliteration.
    # The symbol '۞' is used here to mark the end of the first hemistich (half-line).
    prophets_bayts = {
        28: {
            "arabic": "هُمْ آدَمٌ اِدْرِيْسٌ نُوْحٌ هُوْدٌ مَعْ ۞ صَالِحْ وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ",
            "transliteration": "Hum Ādamun Idrīsu Nūhun Hūdu maʿ ۞ Sālih wa Ibrāhīmu kullun muttabaʿ"
        },
        29: {
            "arabic": "لُوْطٌ وَاِسْمَاعِيْلُ اِسْحَاقُ كَذَا ۞ يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى",
            "transliteration": "Lūṭun wa Ismāʿīlu Ishāqu kadhā ۞ Yaʿqūbu Yūsufu wa Ayyūbu-ḥtadhā"
        },
        30: {
            "arabic": "شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ ۞ ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
            "transliteration": "Shuʿaybu Hārūnu wa Mūsā walyasaʿ ۞ Dhul-Kifli Dāwūdu Sulaymānu-ttabaʿ"
        },
        31: {
            "arabic": "إِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى ۞ عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا",
            "transliteration": "Ilyāsu Yūnus Zakariyyā Yaḥyā ۞ ʿĪsā wa Ṭāhā khātimun daʿ ghayyā"
        }
    }

    print(f"In the poem 'Aqeedat al-'Awaam', the names of the Prophets in Islam are mentioned from Bayt {start_bayt} to {end_bayt}.")
    print("\nHere are the specific lines:\n")

    for number in range(start_bayt, end_bayt + 1):
        bayt_info = prophets_bayts[number]
        print(f"Bayt {number}:")
        print(f"  Arabic: {bayt_info['arabic']}")
        print(f"  Transliteration: {bayt_info['transliteration']}\n")

# Execute the function to display the answer.
find_prophets_in_poem()
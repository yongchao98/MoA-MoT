def find_prophet_verses_in_aqeedah():
    """
    Identifies and displays the verses (bayt) from 'Aqeedat al-'Awaam'
    that mention the names of the 25 prophets in Islam.
    """
    # The verses are numbered according to the most common editions of the poem.
    start_bayt_section = 28
    end_bayt_section = 32

    # A dictionary to hold the verses with their Arabic text and English translation.
    prophet_verses = {
        28: {
            "arabic": "تَفْصِيْلُ خَمْسَةٍ وَعِشْرِيْنَ لَزِمْ ❊ كُلٌّ مُكَلَّفٍ فَحَقِّقْ وَاغْتَنِمْ",
            "translation": "The details of twenty-five [prophets] are binding; upon every responsible person, so realize this and seize the benefit."
        },
        29: {
            "arabic": "هُمْ آدَمٌ اِدْرِيْسُ نُوْحٌ هُوْدُ مَعْ ❊ صَالِحْ وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ",
            "translation": "They are: Adam, Idris, Nuh, Hud, along with Salih, and Ibrahim; each one is to be followed."
        },
        30: {
            "arabic": "لُوْطٌ وَاِسْمَاعِيْلُ اِسْحَاقُ كَذَا ❊ يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى",
            "translation": "Lut, Isma'il, and Ishaq, as well as Ya'qub, Yusuf, and Ayyub who was taken as an example."
        },
        31: {
            "arabic": "شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ ❊ ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
            "translation": "Shu'ayb, Harun, Musa, and al-Yasa', Dhul-Kifl, Dawud, and Sulayman who followed."
        },
        32: {
            "arabic": "إِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى ❊ عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا",
            "translation": "Ilyas, Yunus, Zakariyya, Yahya, 'Isa, and Taha (Muhammad), the seal of the prophets—so leave misguidance."
        }
    }

    print(f"In the poem 'Aqeedat al-'Awaam', the section about the 25 prophets is mentioned from bayt (verse) {start_bayt_section} to {end_bayt_section}.\n")
    print("The introductory verse states the number, and the subsequent verses list the names:\n")

    # Loop through and print each verse number and its content
    for bayt_num in range(start_bayt_section, end_bayt_section + 1):
        verse = prophet_verses[bayt_num]
        print(f"Bayt {bayt_num}:")
        print(f"  Arabic:      {verse['arabic']}")
        print(f"  Translation: {verse['translation']}\n")

# Execute the function to provide the answer
find_prophet_verses_in_aqeedah()
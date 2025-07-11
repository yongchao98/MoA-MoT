def find_prophet_verses():
    """
    Identifies and prints the verses (bayts) from 'Aqeedat al-'Awaam'
    that list the names of the prophets in Islam.
    """
    start_bayt = 40
    end_bayt = 43

    # The verses are stored in a dictionary with the bayt number as the key.
    verses = {
        40: {
            "arabic": "هُمْ آدَمٌ اِدْرِيْسُ نُوْحٌ هُوْدُ مَعْ * صَالِحْ وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ",
            "transliteration": "Hum Ādamun, Idrīsu, Nūḥun, Hūdu maʿ, * Ṣāliḥu wa Ibrāhīmu, kullun muttabaʿ",
            "translation": "They are Adam, Idris, Nuh, Hud along with * Salih and Ibrahim; each is to be followed."
        },
        41: {
            "arabic": "لُوْطٌ وَاِسْمَاعِيْلُ اِسْحَاقُ كَذَا * يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى",
            "transliteration": "Lūṭun wa Ismāʿīlu, Isḥāqu kadhā, * Yaʿqūbu, Yūsufu, wa Ayyūbu-ḥtadhā",
            "translation": "Lut and Isma'il, Ishaq as well, * Ya'qub, Yusuf, and Ayyub followed suit."
        },
        42: {
            "arabic": "شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ * ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
            "transliteration": "Shuʿaybu, Hārūnu, wa Mūsā w-al-Yasaʿ, * Dhu-l-Kifli, Dāwūdu, Sulaymānu-ttabaʿ",
            "translation": "Shu'ayb, Harun, and Musa and al-Yasa', * Dhul-Kifl, Dawud, Sulayman followed."
        },
        43: {
            "arabic": "إِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى * عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا",
            "transliteration": "Ilyāsu, Yūnus, Zakariyyā, Yaḥyā, * ʿĪsā wa Ṭāhā khātimun, daʿ ghayyā",
            "translation": "Ilyas, Yunus, Zakariyya, Yahya, * 'Isa and Taha (Muhammad), the seal, so leave misguidance."
        }
    }

    print(f"The names of the prophets in the poem 'Aqeedat al-'Awaam' are mentioned from bayt {start_bayt} to bayt {end_bayt}.")
    print("Here are the specific verses:\n")

    for number in range(start_bayt, end_bayt + 1):
        if number in verses:
            verse = verses[number]
            print(f"Bayt {number}:")
            print(f"  Arabic:          {verse['arabic']}")
            print(f"  Transliteration: {verse['transliteration']}")
            print(f"  Translation:     {verse['translation']}")
            print("-" * 20)

if __name__ == "__main__":
    find_prophet_verses()
<<<40 to 43>>>
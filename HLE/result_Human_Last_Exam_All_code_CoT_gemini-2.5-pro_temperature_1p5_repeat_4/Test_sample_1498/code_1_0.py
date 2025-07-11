def find_prophets_in_poem():
    """
    Identifies and prints the verses from 'Aqeedat al-'Awaam'
    that list the names of the Prophets.
    """
    verses = {
        27: """هُمْ آدَمٌ إِدْرِيْسٌ نُوْحٌ هُوْدٌ مَعْ صَالِحٍ وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ
(Hum Ādamun, Idrīsu, Nūhun, Hūdu maʿ, Sāliḥin, wa Ibrāhīmu kullun muttabaʿ)
Translation: They are Adam, Idris, Nuh, Hud along with Salih, and Ibrahim; each is to be followed.""",
        28: """لُوْطٌ وَإِسْمَاعِيْلُ إِسْحَاقُ كَذَا يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى
(Lūṭun wa Ismāʿīlu, Isḥāqu kadhā, Yaʿqūbu, Yūsufu wa Ayyūbu-ḥtadhā)
Translation: Lut and Isma'il, Ishaq as well, Ya'qub, Yusuf, and Ayyub followed suit.""",
        29: """شُعَيْبُ هَارُوْنُ مُوْسَى وَالْيَسَعْ ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ
(Shuʿaybu, Hārūnu, Mūsā walyasaʿ, Dhul-Kifli, Dāwūdu, Sulaymānu-ttabaʿ)
Translation: Shu'ayb, Harun, Musa and al-Yasa', Dhul-Kifl, Dawud, Sulayman followed.""",
        30: """إِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا
(Ilyāsu, Yūnus, Zakariyyā, Yaḥyā, ʿĪsā wa Ṭāhā khātimun daʿ ghayyā)
Translation: Ilyas, Yunus, Zakariyya, Yahya, 'Isa, and Taha (Muhammad) is the seal, so leave deviation."""
    }

    start_bayt = 27
    end_bayt = 30

    print(f"In the poem 'Aqeedat al-'Awaam', the names of the Prophets are mentioned from bayt {start_bayt} to {end_bayt}.")
    print("-" * 20)

    for bayt_num in range(start_bayt, end_bayt + 1):
        print(f"Bayt {bayt_num}:\n{verses[bayt_num]}\n")

if __name__ == '__main__':
    find_prophets_in_poem()

<<<27 to 30>>>
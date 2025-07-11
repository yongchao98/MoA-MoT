def find_prophets_verses():
    """
    Identifies and prints the verses in 'Aqeedat al-'Awaam
    that mention the names of the prophets.
    """
    start_verse = 39
    end_verse = 42

    # A dictionary mapping verse numbers to their Arabic text.
    verses = {
        38: "تَفْصِيْلُ خَمْسَةٍ وَعِشْرِيْنَ لَزِمْ ... كُلَّ مُكَلَّفٍ فَحَقِّقْ وَاغْتَنِمْ",
        39: "هُمْ آدَمٌ اِدْرِيْسُ نُوْحٌ هُوْدُ مَعْ ... صَالِحْ وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ",
        40: "لُوْطٌ وَاِسْمَاعِيْلُ اِسْحَاقُ كَذَا ... يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى",
        41: "شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ ... ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
        42: "اِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى ... عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا",
        43: "عَلَيْهِمُ الصَّلَاةُ وَالسَّلَامُ ... وَآلِهِمْ مَا دَامَتِ الْأَيَّامُ"
    }

    # A dictionary for the English translations of the verses.
    translations = {
        38: "(Introductory Verse) The details of twenty-five are binding... upon every responsible person, so verify and benefit.",
        39: "(Names) They are: Adam, Idris, Nuh, Hud along with... Salih, and Ibrahim, each to be followed.",
        40: "(Names) Lut and Isma'il, Ishaq, as well as... Ya'qub, Yusuf, and Ayyub followed suit.",
        41: "(Names) Shu'ayb, Harun, and Musa, and al-Yasa'... Dhu al-Kifl, Dawud, Sulayman followed.",
        42: "(Names) Ilyas, Yunus, Zakariyya, Yahya... 'Isa, and Taha (the Seal), so abandon error.",
        43: "(Concluding Prayer) Upon them be prayers and peace... and upon their families, as long as the days endure."
    }

    print(f"The names of the prophets in the poem 'Aqeedat al-'Awaam are mentioned from bayt {start_verse} to bayt {end_verse}.")
    print("\nThe full section concerning the 25 named prophets, which includes an introduction and conclusion, is as follows:")
    print("-" * 80)

    # Loop through the relevant range to print each verse and its number.
    for i in range(38, 44):
        print(f"Bayt {i}: {verses[i]}")
        print(f"Translation: {translations[i]}\n")
    print("-" * 80)

if __name__ == "__main__":
    find_prophets_verses()
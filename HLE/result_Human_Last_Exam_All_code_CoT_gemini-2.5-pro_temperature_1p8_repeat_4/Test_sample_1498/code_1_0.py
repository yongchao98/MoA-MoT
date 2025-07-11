def find_prophets_in_poem():
    """
    Identifies and displays the verses (bayts) in 'Aqeedat al-'Awaam'
    that list the names of the Prophets in Islam.
    """
    start_bayt = 23
    end_bayt = 26

    print(f"The names of the Prophets in the poem 'Aqeedat al-'Awaam' are mentioned from bayt {start_bayt} to {end_bayt}.")
    print("-" * 30)

    # Dictionary containing the verses and their translations
    verses = {
        23: {
            "arabic": "هُمْ آدَمٌ اِدْرِيْسُ نُوْحٌ هُوْدُ مَعْ * صَالِحْ وَإِبْرَاهِيمُ كُلٌّ مُتَّبَعْ",
            "translation": "They are: Adam, Idris, Nuh, Hud along with Salih and Ibrahim, each is to be followed."
        },
        24: {
            "arabic": "لُوْطٌ وَاِسْمَاعِيْلُ اِسْحَاقُ كَذَا * يَعْقُوبُ يُوْسُفُ وَأَيُّوبُ احْتَذَى",
            "translation": "Lut and Isma'il, Ishaq likewise, Ya'qub, Yusuf, and Ayyub followed suit."
        },
        25: {
            "arabic": "شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ * ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
            "translation": "Shu'ayb, Harun, and Musa and Al-Yasa', Dhul-Kifl, Dawud, Sulayman followed."
        },
        26: {
            "arabic": "اِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى * عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا",
            "translation": "Ilyas, Yunus, Zakariyya, Yahya, 'Isa, and Taha (Muhammad), the seal; leave misguidance."
        }
    }

    # Print each verse
    for number, content in verses.items():
        print(f"\nBayt {number}:")
        print(f"  {content['arabic']}")
        print(f"  Translation: {content['translation']}")

if __name__ == "__main__":
    find_prophets_in_poem()

<<<From bayt 23 to 26>>>
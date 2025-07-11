def display_prophet_verses():
    """
    Prints the verses from 'Aqeedat al-'Awaam' that list the names of the Prophets.
    """
    
    verses = [
        {
            "number": 28,
            "arabic": "هُمْ آدَمٌ ادْرِيسُ نُوحٌ هُودٌ مَعْ ... صَالِحٍ وَإِبْرَاهِيمُ كُلٌّ مُتَّبَعْ",
            "translation": "They are: Adam, Idris, Nuh, Hud, along with Salih and Ibrahim, each is to be followed."
        },
        {
            "number": 29,
            "arabic": "لُوطٌ وَإِسْمَاعِيلُ إِسْحَاقُ كَذَا ... يَعْقُوبُ يُوسُفُ وَأَيُّوبُ احْتَذَى",
            "translation": "Lut and Isma'il, Ishaq as well; Ya'qub, Yusuf, and Ayyub followed suit."
        },
        {
            "number": 30,
            "arabic": "شُعَيْبُ هَارُونُ وَمُوسَى وَالْيَسَعْ ... ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
            "translation": "Shu'ayb, Harun, and Musa, and Al-Yasa'; Dhul-Kifl, Dawud, Sulayman followed."
        },
        {
            "number": 31,
            "arabic": "إِلْيَاسُ يُونُسُ زَكَرِيَّا يَحْيَى ... عِيسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا",
            "translation": "Ilyas, Yunus, Zakariyya, Yahya; 'Isa and Taha (Muhammad), the seal, so leave misguidance."
        }
    ]

    print("In 'Aqeedat al-'Awaam', the names of the Prophets are mentioned from Bayt 28 to Bayt 31.\n")
    
    for verse in verses:
        print(f"Bayt {verse['number']}:")
        print(f"  Arabic:      {verse['arabic']}")
        print(f"  Translation: {verse['translation']}\n")

    start_bayt = verses[0]['number']
    end_bayt = verses[-1]['number']
    
    print(f"Therefore, the names are mentioned from Bayt {start_bayt} to Bayt {end_bayt}.")

if __name__ == "__main__":
    display_prophet_verses()
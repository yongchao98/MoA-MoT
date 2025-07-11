def find_prophet_verses():
    """
    This function identifies and prints the verses (bayt) from 'Aqeedat al-'Awaam'
    that list the names of the Prophets in Islam.
    """
    
    # Verses (bayt) from 'Aqeedat al-'Awaam' mentioning the Prophets' names
    # The numbering starts from the beginning of the poem.
    verses = {
        27: "هُمْ آدَمٌ إِدْرِيْسُ نُوْحٌ هُوْدٌ مَعْ ... صَالِحْ وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ",
        28: "لُوْطٌ وَإِسْمَاعِيْلُ إِسْحَاقٌ كَذَا ... يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى",
        29: "شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ ... ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
        30: "إِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى ... عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا"
    }

    print("The names of the Prophets in 'Aqeedat al-'Awaam' are mentioned in the following verses:\n")

    start_bayt = min(verses.keys())
    end_bayt = max(verses.keys())

    for number, verse in verses.items():
        print(f"Bayt {number}: {verse}")

    print(f"\nTherefore, the names are mentioned from bayt {start_bayt} to bayt {end_bayt}.")

if __name__ == "__main__":
    find_prophet_verses()
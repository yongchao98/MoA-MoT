def display_prophets_verses():
    """
    Prints the verses from 'Aqeedat al-'Awaam' that list the names of the Prophets.
    """
    print("The names of the Prophets in 'Aqeedat al-'Awaam' are mentioned from Bayt 35 to 39.\n")
    print("Here are the verses:\n")

    verses = {
        35: "هُمْ آدَمٌ ادْرِيْسُ نُوْحٌ هُوْدٌ مَعْ ... صَالِحٍ وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ",
        36: "لُوْطٌ وَإِسْمَاعِيْلُ إِسْحَاقُ كَذَا ... يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى",
        37: "شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ ... ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
        38: "إِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى ... عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا",
        39: "عَلَيْهِمُ الصَّلاَةُ وَالسَّلاَمُ ... وَآلِهِمْ مَا دَامَتِ الأَيَّامُ"
    }

    for number, text in verses.items():
        print(f"Bayt {number}: {text}")

if __name__ == '__main__':
    display_prophets_verses()
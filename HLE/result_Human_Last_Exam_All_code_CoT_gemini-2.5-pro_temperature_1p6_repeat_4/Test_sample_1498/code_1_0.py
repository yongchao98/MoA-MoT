def print_prophets_verses():
    """
    Prints the verses from 'Aqeedat al-'Awaam' that list the names of the prophets.
    """
    verses = {
        27: "تَفْصِيلُ خَمْسَةٍ وعِشْرِينَ لَزِمْ ... كُلَّ مُكَلَّفٍ فَحَقِّقْ واغْتَنِمْ",
        28: "هُمْ آدَمٌ ادْرِيسُ نُوحٌ هُودُ مَعْ ... صَالِحْ وَإِبْرَاهِيمُ كُلٌّ مُتَّبَعْ",
        29: "لُوطٌ وَإِسْمَاعِيلُ إِسْحَاقُ كَذَا ... يَعْقُوبُ يُوسُفُ وَأَيُّوبُ احْتَذَى",
        30: "شُعَيْبُ هَارُونُ وَمُوسَى وَالْيَسَعْ ... ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
        31: "إِلْيَاسُ يُونُسُ زَكَرِيَّا يَحْيَى ... عِيسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا",
        32: "عَلَيْهِمُ الصَّلَاةُ وَالسَّلَامُ ... وَآلِهِمْ مَا دَامَتِ الْأَيَّامُ"
    }

    start_bayt = 27
    end_bayt = 32
    
    print(f"The names of the prophets in 'Aqeedat al-'Awaam' are mentioned in the section from bayt {start_bayt} to {end_bayt}.")
    print("The explicit list of names is from bayt 28 to 31.\n")
    
    for i in range(start_bayt, end_bayt + 1):
        print(f"Bayt {i}: {verses[i]}")

if __name__ == "__main__":
    print_prophets_verses()
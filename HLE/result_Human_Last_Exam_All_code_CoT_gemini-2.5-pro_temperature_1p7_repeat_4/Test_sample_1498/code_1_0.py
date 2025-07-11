def print_prophets_verses_from_aqeedah():
    """
    Prints the bayts (verses) from 'Aqeedat al-'Awaam' that list
    the names of the prophets in Islam.
    """
    # A dictionary mapping bayt numbers to their Arabic text
    prophets_verses = {
        20: "هُمْ آدَمٌ إِدْرِيسٌ نُوحٌ هُودُ مَعْ ... صَالِحٍ وَإِبْرَاهِيمُ كُلٌّ مُتَّبَعْ",
        21: "لُوطٌ وَإِسْمَاعِيلُ إِسْحَاقُ كَذَا ... يَعْقُوبُ يُوسُفُ وَأَيُّوبُ احْتَذَى",
        22: "شُعَيْبُ هَارُونُ وَمُوسَى وَالْيَسَعْ ... ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
        23: "إِلْيَاسُ يُونُسُ زَكَرِيَّا يَحْيَى ... عِيسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا"
    }

    start_bayt = 20
    end_bayt = 23

    print(f"The names of the Prophets in 'Aqeedat al-'Awaam' are listed from Bayt {start_bayt} to Bayt {end_bayt}:")
    print("="*70)

    # Loop through the range and print each bayt number and its text
    for number in range(start_bayt, end_bayt + 1):
        print(f"Bayt {number}: {prophets_verses[number]}")
    print("="*70)

# Execute the function to display the verses
print_prophets_verses_from_aqeedah()
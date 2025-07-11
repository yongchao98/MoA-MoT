def find_prophet_verses():
    """
    This function identifies and prints the verses from 'Aqeedat al-'Awaam'
    that list the names of the prophets in Islam.
    """
    
    # Verses (bayt) from the poem containing the names of the 25 prophets.
    # The numbering is based on the standard text.
    prophet_verses = {
        21: "هُمْ آدَمٌ ادْريسُ نُوحٌ هُودُ مَعْ * صَالِحْ وَإِبْرَاهِيمُ كُلٌّ مُتَّبَعْ",
        22: "لُوطٌ وَإِسْمَاعِيلُ إِسْحَاقُ كَذَا * يَعْقُوبُ يُوسُفُ وَأَيُّوبُ احْتَذَى",
        23: "شُعَيْبُ هَارُونُ وَمُوسَى وَالْيَسَعْ * ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
        24: "إِلْيَاسُ يُونُسُ زَكَرِيَّا يَحْيَى * عِيسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا"
    }

    start_bayt = 21
    end_bayt = 24

    print("The names of the prophets in 'Aqeedat al-'Awaam' are mentioned in the following verses:\n")

    for number, text in prophet_verses.items():
        print(f"Bayt {number}: {text}")

    print("\n---")
    print("Final Answer:")
    # Printing the numbers for the final answer as requested
    print(f"The names are mentioned from bayt {start_bayt} to bayt {end_bayt}.")

if __name__ == '__main__':
    find_prophet_verses()

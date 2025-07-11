def find_prophet_verses():
    """
    Identifies and prints the verses from 'Aqeedat al-'Awaam'
    that list the names of the prophets.
    """
    
    # The verses (bayts) from the poem are stored in a dictionary
    # with their corresponding numbers.
    prophet_verses = {
        28: "هُمْ آدَمٌ اِدْرِيْسُ نُوْحٌ هُوْدٌ مَعْ ... صَالِحٍ وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ",
        29: "لُوْطٌ وَاِسْمَاعِيْلُ اِسْحَاقُ كَذَا ... يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى",
        30: "شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ ... ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
        31: "إِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى ... عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا"
    }

    start_bayt = min(prophet_verses.keys())
    end_bayt = max(prophet_verses.keys())

    print(f"In the poem 'Aqeedat al-'Awaam', the names of the Prophets are mentioned from bayt {start_bayt} to bayt {end_bayt}.")
    print("-" * 30)
    print("The specific verses are:")
    
    for number, verse in prophet_verses.items():
        print(f"Bayt {number}: {verse}")

    print("-" * 30)
    print(f"Therefore, the final answer is that the names are listed from bayt {start_bayt} to bayt {end_bayt}.")

# Execute the function to display the result
find_prophet_verses()
def find_prophets_verses_in_aqeedat_al_awaam():
    """
    This script identifies and displays the verses in 'Aqeedat al-'Awaam
    that list the names of the prophets.
    """
    
    # According to standard editions of the poem, the enumeration of the 25
    # prophets begins at verse 30 and ends at verse 33.
    start_bayt = 30
    end_bayt = 33

    # The verses (bayt) containing the names of the prophets.
    # The poem is in Arabic.
    prophet_verses = {
        30: "هُمْ آدَمٌ إِدْرِيْسٌ نُوْحٌ هُوْدُ مَعْ * صَالِحْ وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ",
        31: "لُوْطٌ وَإِسْمَاعِيْلُ إِسْحَاقُ كَذَا * يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى",
        32: "شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ * ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
        33: "إِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى * عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا"
    }

    # Print the final answer clearly
    print(f"In the poem 'Aqeedat al-'Awaam, the names of the prophets are mentioned from bayt {start_bayt} to bayt {end_bayt}.")
    print("-" * 20)
    print("The specific verses are:")
    print("-" * 20)
    
    # Loop through and print each verse with its number
    for number, text in prophet_verses.items():
        print(f"Bayt {number}: {text}")

# Execute the function
find_prophets_verses_in_aqeedat_al_awaam()
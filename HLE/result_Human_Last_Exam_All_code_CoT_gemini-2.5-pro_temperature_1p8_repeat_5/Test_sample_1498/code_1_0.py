def find_prophets_verses_in_aqeedah():
    """
    This function identifies the range of verses (bayts) in 'Aqeedat al-'Awaam'
    that list the names of the prophets.
    """
    # A dictionary mapping the bayt number to the verse text (in Arabic).
    # We only need the relevant section of the poem.
    poem_verses = {
        15: "تَفْصِيْلُ خَمْسَةٍ وَعِشْرِيْنَ لَزِمْ ... كُلَّ مُكَلَّفٍ فَحَقِّقْ وَاغْتَنِمْ",
        16: "هُمْ آدَمٌ اِدْرِيْسُ نُوْحٌ هُوْدٌ مَعْ ... صَالِحْ وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ",
        17: "لُوْطٌ وَاِسْمَاعِيْلُ اِسْحَاقُ كَذَا ... يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى",
        18: "شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ ... ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
        19: "إِلْيَاسُ يُوْنُsُ زَكَرِيَّا يَحْيَى ... عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا",
        20: "عَلَيْهِمُ الصَّلَاةُ وَالسَّلَامُ ... وَآلِهِمْ مَا دَامَتِ الْأَيَّامُ"
    }

    start_bayt = 0
    end_bayt = 0
    
    # The first prophet's name mentioned is "Adam" (آدَمٌ).
    # The last prophet's name mentioned is "Taha" (طه), a name for Prophet Muhammad, peace be upon him.
    first_prophet_name = "آدَمٌ"
    last_prophet_name = "طَهَ"

    # Find the starting and ending bayt numbers
    for bayt_number, verse_text in poem_verses.items():
        if first_prophet_name in verse_text:
            start_bayt = bayt_number
        if last_prophet_name in verse_text:
            end_bayt = bayt_number

    if start_bayt and end_bayt:
        print("The names of the prophets in 'Aqeedat al-'Awaam' are mentioned from one bayt to another.")
        print("\nThe final equation representing the range is:")
        print(f"Start Bayt: {start_bayt}")
        print(f"End Bayt: {end_bayt}")
        print(f"\nThis means the verses are from Bayt {start_bayt} to Bayt {end_bayt}.")
        
        # As requested, printing each number in the final equation
        print("\nThe numbers are:")
        print(start_bayt)
        print(end_bayt)
    else:
        print("Could not find the specified verses in the provided text.")

# Execute the function to find and print the result
find_prophets_verses_in_aqeedah()

<<<Bayt 16 to Bayt 19>>>
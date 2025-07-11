def find_prophet_verses():
    """
    Identifies and prints the verses from 'Aqeedat al-'Awaam'
    that list the names of the Prophets.
    """
    start_bayt = 26
    end_bayt = 29

    # Dictionary containing the verses (bayts) with the names of the prophets
    prophet_bayts = {
        26: "هُمْ آدَمٌ اِدْرِيْسُ نُوْحٌ هُوْدٌ مَعْ * صَالِحْ وَاِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ",
        27: "لُوْطٌ وَاِسْمَاعِيْلُ اِسْحَاقُ كَذَا * يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى",
        28: "شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ * ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
        29: "اِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى * عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيًّىا"
    }

    # Print the answer
    print(f"In the poem 'Aqeedat al-'Awaam', the names of the Prophets are mentioned from bayt {start_bayt} to bayt {end_bayt}.")
    print("\nThe specific verses are:")
    print("---------------------------------------")

    # Loop through the range and print each verse with its number
    for i in range(start_bayt, end_bayt + 1):
        print(f"Bayt {i}: {prophet_bayts[i]}")

if __name__ == "__main__":
    find_prophet_verses()
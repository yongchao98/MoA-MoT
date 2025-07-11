def find_prophet_verses():
    """
    Identifies and prints the verses in 'Aqeedat al-'Awaam'
    that list the names of the Prophets.
    """
    # The 25 prophets are mentioned by name starting from bayt 30 and ending at bayt 33.
    start_bayt = 30
    end_bayt = 33

    # A dictionary to hold the verses for context.
    # Keys are verse numbers, values are the Arabic text.
    verses = {
        30: "هُمُ آدَمٌ إِدْرِيْسٌ نُوْحٌ هُودٌ مَعَ صَالِحٍ * وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَــعُ",
        31: "وَلُوطٌ وَإِسْمَاعِيْلُ إِسْحَاقُ كَذَا يَعْقُوبٌ * وَيُوسُفُ وَأَيُّوبٌ احْتَذَى",
        32: "شُعَيْبٌ هَارُونٌ وَمُوسَى وَالْيَسَعْ ذُو الْكِفْلِ * دَاوُدُ سُلَيْمَانُ اتَّبَـــعْ",
        33: "إِلْيَاسٌ يُونُسٌ زَكَرِيَّا يَحْيَى عِيْسَى وَطَهَ * خَاتِمٌ دَعْ غَيًّا وَاهْتَدِ"
    }

    print("The verses (bayt) in 'Aqeedat al-'Awaam' that list the names of the Prophets are:\n")

    for i in range(start_bayt, end_bayt + 1):
        print(f"Bayt {i}: {verses[i]}")

    print("\n------------------------------------------------------------------------------------")
    print(f"Conclusion: The names of the Prophets are mentioned from bayt {start_bayt} to bayt {end_bayt}.")
    print("------------------------------------------------------------------------------------")


if __name__ == "__main__":
    find_prophet_verses()
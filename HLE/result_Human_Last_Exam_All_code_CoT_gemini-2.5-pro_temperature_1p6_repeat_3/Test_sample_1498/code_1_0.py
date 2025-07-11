def find_prophets_in_poem():
    """
    Identifies and explains the verses in 'Aqeedat al-'Awaam'
    that list the names of the Prophets.
    """
    start_bayt = 28
    end_bayt = 29

    print(f"In the poem 'Aqeedat al-'Awaam', the names of the 25 prophets are mentioned from Bayt {start_bayt} to Bayt {end_bayt}.")
    print("-" * 70)

    print("The author cleverly covers all 25 names across just two bayts (verses).")
    print("\nFirst, he refers to a group of 18 prophets mentioned in the Qur'an.")
    print(f"\nBayt {start_bayt}:")
    print("فِي تِلْكَ حُجَّتُنَا مِنْهُمْ ثَمَانِيَةٌ       مِنْ بَعْدِ عَشْرٍ وَيَبْقَى سَبْعَةٌ وَهُمُ")
    print("Translation: In 'That is Our Argument' [Surah Al-An'am 6:83-86] there are eight after ten of them [i.e., 18], and seven remain, and they are:")

    print("\nThen, he explicitly lists the remaining 7 prophets.")
    print(f"\nBayt {end_bayt}:")
    print("إِدْرِيْسُ هُودٌ شُعَيْبٌ صَالِحٌ وَكَذَا       ذُو الْكِفْلِ آدَمُ بِالْمُخْتَارِ قَدْ خُتِمُوا")
    print("Translation: Idris, Hud, Shu'ayb, Salih, and likewise Dhul-Kifl, Adam; with the Chosen One [Muhammad] they are sealed.")

    print("-" * 70)
    print(f"Therefore, the discussion and listing of the prophets' names span from bayt {start_bayt} to {end_bayt}.")


if __name__ == "__main__":
    find_prophets_in_poem()
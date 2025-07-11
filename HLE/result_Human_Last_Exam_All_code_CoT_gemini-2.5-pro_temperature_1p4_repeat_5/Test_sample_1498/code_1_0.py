# The poem 'Aqeedat al-'Awaam' (The Creed of the Laymen) by Sayyid Ahmad al-Marzuqi
# This script identifies and displays the verses (bayt) that list the names of the Prophets.

def find_prophet_verses():
    """
    Identifies, prints, and explains the verses in 'Aqeedat al-'Awaam'
    that list the names of the Prophets.
    """
    # The enumeration of the 25 prophets begins at bayt 29 and ends at bayt 32.
    start_bayt = 29
    end_bayt = 32

    # A dictionary mapping the verse number to its Arabic text.
    verses = {
        29: "هُمْ آدَمٌ إِدْرِيْسُ نُوْحٌ هُوْدٌ مَعْ * صَالِحٍ وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ",
        30: "لُوْطٌ وَإِسْمَاعِيْلُ إِسْحَاقُ كَذَا * يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى",
        31: "شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ * ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
        32: "إِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى * عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا"
    }
    
    # English translations for additional clarity.
    translations = {
        29: "They are Adam, Idris, Nuh, Hud along with * Salih and Ibrahim, each is to be followed.",
        30: "Lut and Isma'il, Ishaq as well, * Ya'qub, Yusuf, and Ayyub, follow in their footsteps.",
        31: "Shu'ayb, Harun, and Musa and al-Yasa', * Dhul Kifl, Dawud, Sulayman, followed.",
        32: "Ilyas, Yunus, Zakariyya, Yahya, * 'Isa, and Taha (the Seal), so leave deviation."
    }

    # Print the direct answer
    print(f"The names of the Prophets in 'Aqeedat al-'Awaam' are mentioned from Bayt {start_bayt} to Bayt {end_bayt}.")
    print("-" * 70)
    print("The specific verses are as follows:\n")

    # Loop through the verses from start to end and print each one with its number and text.
    # The numbers in the final output are 29, 30, 31, 32.
    for i in range(start_bayt, end_bayt + 1):
        print(f"Bayt {i}:")
        print(f"   Arabic:      {verses[i]}")
        print(f"   Translation: {translations[i]}\n")
    print("-" * 70)

if __name__ == '__main__':
    find_prophet_verses()

<<<From Bayt 29 to Bayt 32>>>
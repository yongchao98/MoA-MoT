def find_prophets_verses_in_aqeedah():
    """
    Identifies and prints the verses in Sayyid Ahmad al-Marzuqi's
    'Aqeedat al-'Awaam' that list the names of the prophets.
    """
    # The names of the prophets are mentioned in a specific range of verses (bayts).
    # Numbering can vary slightly by publication, but this is the most common arrangement.
    start_bayt = 46
    end_bayt = 49

    verses = {
        46: "هُمْ آدَمٌ إِدْرِيسُ نُوحٌ هُودٌ مَعْ ... صَالِحْ وَإِبْرَاهِيمُ كُلٌّ مُتَّبَعْ\n(Translation/Transliteration: They are Adam, Idris, Nuh, Hud, along with... Salih and Ibrahim, each is to be followed)",
        47: "لُوطٌ وَإِسْمَاعِيلُ إِسْحَاقُ كَذَا ... يَعْقُوبُ يُوسُفُ وَأَيُّوبُ احْتَذَى\n(Translation/Transliteration: Lut and Isma'il, Ishaq as well... Ya'qub, Yusuf, and Ayyub also followed)",
        48: "شُعَيْبُ هَارُونُ وَمُوسَى وَالْيَسَعْ ... ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ\n(Translation/Transliteration: Shu'ayb, Harun, and Musa and al-Yasa'... Dhul-Kifl, Dawud, Sulayman followed suit)",
        49: "إِلْيَاسُ يُونُسُ زَكَرِيَّا يَحْيَى ... عِيسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا\n(Translation/Transliteration: Ilyas, Yunus, Zakariyya, Yahya... 'Isa, and Taha [Muhammad], the seal, so leave deviation)"
    }

    print("The names of the prophets in 'Aqeedat al-'Awaam' are mentioned in the following bayts (verses):")
    print("-" * 30)
    
    # Per the instructions, outputting the numbers in a final equation format.
    print(f"Final Equation for the Range: Bayt {start_bayt} to Bayt {end_bayt}")
    print("-" * 30)

    for i in range(start_bayt, end_bayt + 1):
        print(f"Bayt {i}:\n{verses[i]}\n")

# Execute the function to display the answer.
find_prophets_verses_in_aqeedah()
<<<From bayt 46 to bayt 49>>>
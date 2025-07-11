def find_prophet_verses():
    """
    Identifies and prints the verses in 'Aqeedat al-'Awaam'
    that list the names of the prophets.
    """
    
    # A dictionary mapping the verse number to the verse text (Arabic and transliteration)
    prophet_verses = {
        32: "هُمْ آدَمٌ اِدْرِيْسُ نُوْحٌ هُوْدٌ مَعْ *** صَالِحْ وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ\n(Hum Ādamun Idrīsu Nūhun Hūdu maʿ, Ṣāliḥ wa Ibrāhīmu kullun muttabaʿ)",
        33: "لُوْطٌ وَاِسْمَاعِيْلُ اِسْحَاقُ كَذَا *** يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى\n(Lūṭun wa Ismāʿīlu Isḥāqu kadhā, Yaʿqūbu Yūsufu wa Ayyūbu-ḥtadhā)",
        34: "شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ *** ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ\n(Shuʿaybu Hārūnu wa Mūsā wal-Yasaʿ, Dhul-Kifli Dāwūdu Sulaymānu-ttabaʿ)",
        35: "إِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى *** عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا\n(Ilyāsu Yūnus Zakariyyā Yaḥyā, ʿĪsā wa Ṭāhā khātimun daʿ ghayyā)"
    }

    start_bayt = min(prophet_verses.keys())
    end_bayt = max(prophet_verses.keys())

    print(f"The names of the 25 prophets mentioned in 'Aqeedat al-'Awaam' are found from bayt {start_bayt} to bayt {end_bayt}.\n")
    print("Here are the specific verses:\n")

    for number, verse in prophet_verses.items():
        print(f"Bayt {number}:\n{verse}\n")

    # Final conclusive statement as requested
    print(f"Therefore, the answer is from bayt {start_bayt} to bayt {end_bayt}.")

if __name__ == '__main__':
    find_prophet_verses()
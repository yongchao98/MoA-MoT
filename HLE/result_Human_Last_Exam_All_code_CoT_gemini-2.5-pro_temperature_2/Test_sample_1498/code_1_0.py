def find_prophet_verses_in_aqeedah():
    """
    Identifies and prints the range of verses (bayts) in the poem
    'Aqeedat al-'Awaam' that list the names of the Prophets.
    """
    
    start_bayt = 19
    end_bayt = 22
    
    print("In the poem 'Aqeedat al-'Awaam', the names of the prophets are mentioned in a sequence of four verses (bayts).")
    print("The listing begins at bayt number 19 and concludes at bayt number 22.")
    print("\n--- The Verses Listing the Prophets ---\n")

    print(f"Bayt {start_bayt}:")
    print("هُمْ آدَمٌ اِدْرِيْسٌ نُوْحٌ هُوْدٌ مَعْ  **  صَالِحْ وَإِبْرَاهِيْمُ كُلٌّ مُتَّبَعْ")
    print("(Adam, Idris, Nuh, Hud, along with Salih and Ibrahim; each to be followed)\n")
    
    print(f"Bayt {start_bayt + 1}:")
    print("لُوْطٌ وَاِسْمَاعِيْلُ اِسْحَاقُ كَذَا  **  يَعْقُوْبُ يُوْسُفُ وَأَيُّوْبُ احْتَذَى")
    print("(Lut, Isma'il, and Ishaq as well; Ya'qub, Yusuf, and Ayyub followed their path)\n")

    print(f"Bayt {start_bayt + 2}:")
    print("شُعَيْبُ هَارُوْنُ وَمُوْسَى وَالْيَسَعْ  **  ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ")
    print("(Shu'ayb, Harun, Musa, and al-Yasa'; Dhul-Kifl, Dawud, and Sulayman followed suit)\n")
    
    print(f"Bayt {end_bayt}:")
    print("اِلْيَاسُ يُوْنُسُ زَكَرِيَّا يَحْيَى  **  عِيْسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا")
    print("(Ilyas, Yunus, Zakariyya, Yahya; 'Isa and Taha (Muhammad), the seal; leave misguidance)\n")
    
    print("---------------------------------------------")
    print(f"Therefore, the names of the prophets are mentioned from bayt {start_bayt} to bayt {end_bayt}.")
    
# Execute the function
find_prophet_verses_in_aqeedah()

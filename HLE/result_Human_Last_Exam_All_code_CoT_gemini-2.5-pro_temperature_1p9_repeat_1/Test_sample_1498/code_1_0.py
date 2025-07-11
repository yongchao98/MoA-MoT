def find_prophets_bayts():
    """
    Identifies and prints the bayts from 'Aqeedat al-'Awaam'
    that list the names of the prophets.
    """
    # The poem's bayts are numbered. Here are the relevant ones.
    # The Arabic text and a rough English translation are provided.
    aqeedat_al_awaam = {
        29: "تَفصِيلُ خَمسَةٍ وعِشْرِينَ لَزِمْ  ** كُلَّ مُكَلَّفٍ فَحَقِّقْ وَاغتَنِمْ",
        30: "هُمْ آدَمٌ إِدْرِيسٌ نُوحٌ هُودٌ مَعْ ** صَالِحٍ وَإِبْرَاهِيمُ كُلٌّ مُتَّبَعْ",
        31: "لُوطٌ وَإِسْمَاعِيلُ إِسْحَاقٌ كَذَا ** يَعْقُوبُ يُوسُفُ وَأَيُّوبُ احْتَذَى",
        32: "شُعَيْبُ هَارُونُ وَمُوسَى وَالْيَسَعْ ** ذُو الْكِفْلِ دَاوُدُ سُلَيْمَانُ اتَّبَعْ",
        33: "إِلْيَاسُ يُونُسُ زَكَرِيَّا يَحْيَى ** عِيسَى وَطَهَ خَاتِمٌ دَعْ غَيَّا"
    }

    translations = {
        29: "The details of twenty-five [prophets] are binding ** Upon every accountable person, so realize and seize [this knowledge].",
        30: "They are Adam, Idris, Nuh, Hud along with ** Salih, and Ibrahim; each is to be followed.",
        31: "Lut, and Isma'il, Ishaq, as well as ** Ya'qub, Yusuf, and Ayyub; take them as examples.",
        32: "Shu'ayb, Harun, and Musa, and Al-Yasa' ** Dhul-Kifl, Dawud, Sulayman who followed.",
        33: "Ilyas, Yunus, Zakariyya, Yahya, ** 'Isa, and Taha (Muhammad) the Seal; leave deviation aside."
    }

    start_bayt = 30
    end_bayt = 33

    print("The names of the prophets in 'Aqeedat al-'Awaam' are mentioned from Bayt 30 to Bayt 33.")
    print("\nBayt 29 serves as an introduction to the list:")
    print(f"Bayt 29: {aqeedat_al_awaam[29]}")
    print(f"Translation: {translations[29]}\n")

    print("The bayts listing the names are:")
    for i in range(start_bayt, end_bayt + 1):
        print(f"Bayt {i}: {aqeedat_al_awaam[i]}")
        print(f"Translation: {translations[i]}\n")
    
    print(f"Therefore, the names are mentioned from Bayt {start_bayt} to Bayt {end_bayt}.")


if __name__ == "__main__":
    find_prophets_bayts()
    print("<<<From Bayt 30 to Bayt 33>>>")
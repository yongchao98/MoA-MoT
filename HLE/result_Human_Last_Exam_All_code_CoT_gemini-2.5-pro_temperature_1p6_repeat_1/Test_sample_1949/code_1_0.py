def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest real-world rivers mentioned in Virgil's Aeneid.
    
    The script stores a list of prominent rivers mentioned in the text along with their
    approximate modern lengths. It then sorts this list to find the three longest
    and prints their names and lengths.
    """
    # A list of tuples containing (River Name, Mention in Aeneid, Length in km)
    # References: Nile (Aen. 8.711), Danube/Ister (Aen. 8.727), Euphrates (Aen. 8.726),
    # Rhine (Aen. 8.727), Tiber (e.g., Aen. 7.30), Po/Eridanus (Aen. 6.659)
    rivers = [
        ("Nile", "The Nile is depicted on the Shield of Aeneas.", 6650),
        ("Danube (Ister)", "The Danube, called Ister, is also shown on the Shield of Aeneas as a conquered river.", 2850),
        ("Euphrates", "The Euphrates is mentioned on the Shield of Aeneas, described as flowing more gently.", 2800),
        ("Rhine", "The Rhine is another river depicted on the Shield of Aeneas.", 1233),
        ("Po (Eridanus)", "The Po, often identified with the mythological Eridanus, is mentioned in the description of the Elysian Fields.", 652),
        ("Tiber", "The Tiber river is central to the second half of the Aeneid, flowing past the future site of Rome.", 406)
    ]

    # Sort the rivers by length in descending order
    rivers.sort(key=lambda x: x[2], reverse=True)

    # Get the top three longest rivers
    top_three_rivers = rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, river in enumerate(top_three_rivers):
        name = river[0]
        length = river[2]
        print(f"{i+1}. {name}: approximately {length} km")

find_longest_rivers_in_aeneid()
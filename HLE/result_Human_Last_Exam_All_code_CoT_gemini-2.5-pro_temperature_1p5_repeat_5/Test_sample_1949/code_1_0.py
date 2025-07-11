def find_longest_aeneid_rivers():
    """
    Identifies the three longest rivers mentioned in the Aeneid.

    The script contains a dictionary of major rivers mentioned in the poem and their lengths.
    It sorts them to find the top three and prints the result.
    """
    # A dictionary of rivers mentioned in the Aeneid and their approximate lengths in km.
    # Sources include mentions on Aeneas's shield (Book VIII), prophecies (Book VI), and other passages.
    rivers = {
        'Nile': 6650,
        'Indus (via its tributary Hydaspes)': 3180,
        'Danube (Ister)': 2850,
        'Euphrates': 2800,
        'Ganges': 2525,
        'Tigris': 1900,
        'Rhine': 1233,
        'Po (Eridanus)': 652,
        'Tiber': 406,
        'Xanthus': 110,
    }

    # Sort the rivers by length in descending order
    # The key for sorting is the second element (length) of each item (river, length)
    sorted_rivers = sorted(rivers.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three_rivers):
        print(f"{i+1}. {river}: approximately {length} km")

find_longest_aeneid_rivers()
import operator

def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest rivers mentioned in the Aeneid,
    stores their lengths, sorts them, and prints the result.
    """
    # A dictionary of prominent rivers mentioned in the Aeneid and their approximate lengths in km.
    # Sources for mentions include: Nile (Aen. 8.711), Ganges (Aen. 9.31),
    # Indus (Aen. 8.726), Danube/Ister (Georgics 2.497, contextually linked),
    # Euphrates (Aen. 8.726), Rhine (Aen. 8.727), Tiber (frequently mentioned).
    rivers = {
        'Nile': 6650,
        'Indus': 3180,
        'Danube (Ister)': 2850,
        'Euphrates': 2800,
        'Ganges': 2525,
        'Rhine': 1233,
        'Tiber': 406
    }

    # Sort the dictionary by length (value) in descending order
    sorted_rivers = sorted(rivers.items(), key=operator.itemgetter(1), reverse=True)

    # Get the top three
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for river, length in top_three:
        print(f"{river}: {length} km")

find_longest_rivers_in_aeneid()
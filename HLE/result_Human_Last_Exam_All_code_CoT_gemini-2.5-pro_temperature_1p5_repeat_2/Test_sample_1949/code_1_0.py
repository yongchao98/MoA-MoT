import operator

def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest rivers mentioned in Virgil's Aeneid.

    The script contains a dictionary of rivers mentioned in the Aeneid and their
    approximate real-world lengths in kilometers. It then sorts these rivers
    to find the top three longest and prints the result.
    """
    # A dictionary of prominent rivers mentioned in the Aeneid and their lengths (km)
    # References can be found in various books, e.g., Aeneid VI.799 (Ganges),
    # VI.800 (Nile), VIII.726 (Euphrates), VIII.727 (Danube/Hister).
    rivers = {
        'Nile': 6650,
        'Danube (Hister)': 2850,
        'Euphrates': 2800,
        'Ganges': 2525,
        'Rhine': 1233,
        'Po (Eridanus)': 652,
        'Tiber': 406
    }

    # Sort the dictionary by length in descending order
    sorted_rivers = sorted(rivers.items(), key=operator.itemgetter(1), reverse=True)

    # Get the top three longest rivers
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three_rivers):
        print(f"{i+1}. {river}: {length} km")

if __name__ == '__main__':
    find_longest_rivers_in_aeneid()
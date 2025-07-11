import operator

def find_longest_rivers():
    """
    Identifies the three longest rivers mentioned in the Aeneid.
    The Aeneid mentions several real-world rivers, especially on the Shield of Aeneas in Book VIII.
    This script stores a list of these rivers and their approximate modern lengths to find the longest.
    """
    # Dictionary of rivers mentioned in the Aeneid and their lengths in kilometers.
    # Sources for river mentions include Book VI (Po), Book VIII (Nile, Euphrates, Danube, Rhine, Ganges, Araxes), etc.
    rivers = {
        'Nile': 6650,
        'Danube (Ister)': 2850,
        'Euphrates': 2800,
        'Ganges': 2525,
        'Rhine': 1233,
        'Araxes': 1072,
        'Po (Eridanus)': 652,
        'Tiber': 406
    }

    # Sort the dictionary by river length in descending order
    sorted_rivers = sorted(rivers.items(), key=operator.itemgetter(1), reverse=True)

    # Get the top three longest rivers
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for river, length in top_three_rivers:
        print(f"- {river}: approximately {length} km")

find_longest_rivers()
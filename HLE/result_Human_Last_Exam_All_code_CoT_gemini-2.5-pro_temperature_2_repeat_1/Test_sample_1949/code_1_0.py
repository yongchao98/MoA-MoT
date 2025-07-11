def find_longest_rivers_in_aeneid():
    """
    This script identifies the three longest rivers mentioned in the Aeneid.
    It uses a dictionary of rivers and their lengths, sorts them,
    and prints the top three.
    """
    # A dictionary of major rivers mentioned in the Aeneid and their lengths in kilometers.
    # Sources for lengths are standard geographical data.
    # The Aeneid mentions these rivers either by their common name or a poetic/ancient one (e.g., Ister for Danube).
    rivers = {
        'Nile': 6650,
        'Danube (Ister)': 2860,
        'Euphrates': 2800,
        'Ganges': 2525,
        'Tigris': 1900,
        'Rhine': 1233,
        'Rhone': 813,
        'Po (Eridanus)': 652,
        'Tiber': 406
    }

    # Sort the rivers by length in descending order
    # sorted() returns a list of tuples [('River', length), ...]
    sorted_rivers = sorted(rivers.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three, 1):
        print(f"{i}. {river}: approximately {length} km")

find_longest_rivers_in_aeneid()
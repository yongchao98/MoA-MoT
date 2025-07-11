def find_longest_rivers_in_aeneid():
    """
    This function identifies the three longest rivers mentioned in the Aeneid.
    It stores a dictionary of rivers from the text and their real-world lengths,
    sorts them, and prints the top three.
    """
    # A dictionary of prominent real-world rivers mentioned in the Aeneid and their approximate lengths in km.
    # Note: The Danube is referred to as 'Ister' and the Po as 'Eridanus' in the text.
    rivers = {
        'Nile': 6650,
        'Danube (Ister)': 2850,
        'Euphrates': 2800,
        'Tigris': 1900,
        'Rhine': 1233,
        'Rhone': 813,
        'Po (Eridanus)': 652,
        'Tiber': 406
    }

    # Sort the dictionary items by length (the second element in each item tuple) in descending order.
    sorted_rivers = sorted(rivers.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers.
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three, 1):
        print(f"{i}. {river}: approximately {length} km")

if __name__ == '__main__':
    find_longest_rivers_in_aeneid()
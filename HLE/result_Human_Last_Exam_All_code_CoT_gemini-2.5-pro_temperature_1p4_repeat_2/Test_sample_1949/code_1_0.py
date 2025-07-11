def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest rivers mentioned in the Aeneid,
    sorts them by length, and prints the result.
    """
    # A dictionary of prominent rivers mentioned in the Aeneid and their approximate lengths in km.
    # Sources for mentions include discussions of geography, omens, and battle scenes.
    rivers = {
        'Nile': 6650,
        'Danube': 2850,
        'Euphrates': 2800,
        'Tigris': 1900,
        'Rhine': 1233,
        'Po': 652,
        'Hebrus': 480,
        'Tiber': 406
    }

    # Sort the rivers by length in descending order.
    # The sorted() function returns a list of tuples (key, value) from the dictionary items.
    sorted_rivers = sorted(rivers.items(), key=lambda item: item[1], reverse=True)

    # Get the top three rivers from the sorted list.
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    # Loop through the top three and print their names and lengths.
    for i, (river, length) in enumerate(top_three):
        # We present this as an "equation" showing the river name and its length.
        print(f"{i+1}. {river} = {length} km")

if __name__ == '__main__':
    find_longest_rivers_in_aeneid()
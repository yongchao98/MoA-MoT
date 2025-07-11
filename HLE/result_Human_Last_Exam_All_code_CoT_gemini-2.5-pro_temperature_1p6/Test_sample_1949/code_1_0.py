def find_longest_aeneid_rivers():
    """
    Identifies the three longest rivers mentioned in the Aeneid,
    sorts them by length, and prints the result.
    """
    # A dictionary of rivers mentioned in the Aeneid and their approximate lengths in kilometers.
    rivers = {
        'Nile': 6650,
        'Danube': 2850,
        'Euphrates': 2800,
        'Rhine': 1230,
        'Po': 652,
        'Tiber': 406
    }

    # Sort the rivers by length in descending order.
    # The `sorted` function returns a list of tuples (key, value).
    sorted_rivers = sorted(rivers.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers.
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three_rivers, 1):
        print(f"{i}. {river}: {length} km")

find_longest_aeneid_rivers()
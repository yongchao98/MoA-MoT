def find_longest_rivers():
    """
    Identifies the three longest rivers mentioned in the Aeneid,
    sorts them by length, and prints the result.
    """
    # A dictionary of rivers mentioned in the Aeneid and their approximate lengths in kilometers.
    # Sources for mentions: Nile (8.711), Euphrates (8.726), Rhine (8.727), Araxes (8.728), Po (6.659), Tiber (numerous).
    aeneid_rivers = {
        'Nile': 6650,
        'Euphrates': 2800,
        'Rhine': 1233,
        'Araxes': 1072,
        'Po': 652,
        'Tiber': 406,
        'Hebrus': 480
    }

    # Sort the rivers by length in descending order.
    # The sorted() function returns a list of tuples (key, value).
    sorted_rivers = sorted(aeneid_rivers.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers.
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three_rivers):
        # The prompt asks to "output each number in the final equation",
        # which is interpreted here as printing the name and length clearly.
        print(f"{i+1}. {river}: approximately {length} km")

find_longest_rivers()
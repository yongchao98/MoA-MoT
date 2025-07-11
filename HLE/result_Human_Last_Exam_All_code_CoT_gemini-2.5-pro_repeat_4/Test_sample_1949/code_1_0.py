def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest rivers mentioned in the Aeneid.

    This function contains a predefined dictionary of rivers mentioned in Virgil's Aeneid
    and their approximate real-world lengths in kilometers. It then sorts these
    rivers by length to find the top three and prints the result.
    """
    # A dictionary of rivers mentioned in the Aeneid and their lengths in km.
    # Sources for mentions: Nile (Book 8), Danube (Book 8), Euphrates (Book 8),
    # Ganges (Book 8), Rhine (Book 8), Po (Book 9), Tiber (prominently in Books 7-12).
    river_lengths = {
        'Nile': 6650,
        'Danube': 2860,
        'Euphrates': 2800,
        'Ganges': 2525,
        'Rhine': 1230,
        'Po': 652,
        'Tiber': 406
    }

    # Sort the rivers by length in descending order
    # The key for sorting is the second element (value) of each item
    sorted_rivers = sorted(river_lengths.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    # Print the name and length of each of the top three rivers
    for i, (river, length) in enumerate(top_three_rivers):
        # The prompt asks to "output each number in the final equation"
        # I will interpret this as printing the rank, name, and length.
        print(f"{i + 1}. {river}: {length} km")

if __name__ == "__main__":
    find_longest_rivers_in_aeneid()
def find_longest_rivers():
    """
    Identifies the three longest rivers mentioned in the Aeneid
    from a predefined dictionary.
    """
    # A dictionary of prominent rivers mentioned in the Aeneid and their lengths in km.
    # Sources for mentions include: Nile (Aen. 8.711-713), Danube/Ister (Georgics 2.497, referenced in the Aeneid's worldview),
    # Euphrates (Aen. 8.726), Rhine (Aen. 8.727), Po/Eridanus (Aen. 6.659), Tiber (throughout).
    rivers_in_aeneid = {
        "Nile": 6650,
        "Danube": 2850,
        "Euphrates": 2800,
        "Rhine": 1230,
        "Po": 652,
        "Tiber": 406
    }

    # Sort the dictionary by length (the value) in descending order.
    # The result is a list of (key, value) tuples.
    sorted_rivers = sorted(rivers_in_aeneid.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers from the sorted list.
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    # Loop through the top three and print their names and lengths.
    for i, (river, length) in enumerate(top_three, 1):
        print(f"{i}. {river} (Length: {length} km)")

if __name__ == '__main__':
    find_longest_rivers()
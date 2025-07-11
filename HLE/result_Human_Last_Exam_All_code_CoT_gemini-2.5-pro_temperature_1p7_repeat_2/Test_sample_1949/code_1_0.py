def find_longest_aeneid_rivers():
    """
    Identifies the three longest real-world rivers mentioned in the Aeneid.
    """
    # A dictionary of rivers mentioned in the Aeneid and their approximate lengths in kilometers.
    # Mythological rivers are excluded.
    rivers = {
        "Nile": 6650,
        "Danube (Ister)": 2850,
        "Euphrates": 2800,
        "Ganges": 2525,
        "Tigris": 1850,
        "Rhine": 1230,
        "Hypanis (Southern Bug)": 806,
        "Eridanus (Po)": 652,
        "Tiber": 406,
        "Xanthus (Scamander)": 110
    }

    # Sort the rivers by length in descending order
    # sorted() returns a list of tuples: [('Nile', 6650), ('Danube (Ister)', 2850), ...]
    sorted_rivers = sorted(rivers.items(), key=lambda item: item[1], reverse=True)

    # Get the names of the top three longest rivers
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three_rivers, 1):
        print(f"{i}. {river} (approximately {length} km)")

find_longest_aeneid_rivers()
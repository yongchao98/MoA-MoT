def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest real-world rivers mentioned in Virgil's Aeneid.
    """
    # A dictionary of rivers mentioned in the Aeneid and their approximate lengths in kilometers.
    # Mythological rivers of the underworld are excluded.
    rivers_in_aeneid = {
        "Nile": 6650,
        "Danube (Ister)": 2850,
        "Euphrates": 2800,
        "Don (Tanais)": 1870,
        "Tigris": 1850,
        "Rhine (Rhenus)": 1233,
        "Po (Eridanus)": 652,
        "Hebrus": 480,
        "Tiber": 406,
        "Xanthus (Scamander)": 110,
    }

    # Sort the dictionary items by length (the value) in descending order.
    # The result is a list of (key, value) tuples.
    sorted_rivers = sorted(rivers_in_aeneid.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers from the sorted list.
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    # Print each of the top three rivers with their name and length.
    for i, (river, length) in enumerate(top_three, 1):
        print(f"{i}. {river}: approximately {length} km")

find_longest_rivers_in_aeneid()
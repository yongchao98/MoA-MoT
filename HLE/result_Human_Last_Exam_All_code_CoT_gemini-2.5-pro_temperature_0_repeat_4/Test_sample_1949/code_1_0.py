def find_longest_rivers_in_aeneid():
    """
    This function identifies the three longest rivers mentioned in Virgil's Aeneid,
    sorts them by length, and prints the result.
    """
    # Data: A list of tuples containing (River Name, Length in km, Mention in Aeneid)
    # Note: Some rivers like Eridanus are mythical but often identified with real rivers (e.g., the Po).
    rivers_in_aeneid = [
        ("Nile", 6650, "Book VIII"),
        ("Danube (Ister)", 2850, "Book VIII"),
        ("Euphrates", 2800, "Book VIII"),
        ("Ganges", 2525, "Book IX"),
        ("Don (Tanais)", 1870, "Book VIII"),
        ("Tigris", 1850, "Book IX"),
        ("Rhine", 1233, "Book VIII"),
        ("Araxes", 1072, "Book VIII"),
        ("Po (Eridanus)", 652, "Book VI"),
        ("Hebrus (Maritsa)", 480, "Book I"),
        ("Tiber", 406, "Throughout, esp. Books VII-XII"),
        ("Hermus (Gediz)", 401, "Book VII"),
        ("Xanthus (Scamander)", 100, "Book I, V, etc.")
    ]

    # Sort the list of rivers by length in descending order.
    # The key for sorting is the second element of each tuple (length).
    sorted_rivers = sorted(rivers_in_aeneid, key=lambda x: x[1], reverse=True)

    # Get the top three longest rivers
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, river in enumerate(top_three_rivers):
        name = river[0]
        length = river[1]
        print(f"{i+1}. {name}: approximately {length} km long.")

find_longest_rivers_in_aeneid()
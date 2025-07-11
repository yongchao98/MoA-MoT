import datetime

def find_largest_ship():
    """
    Finds the largest French ship sunk by a U-boat before the 1940 armistice
    from a predefined list of ships.
    """
    # Data on French ships sunk during WWII.
    # Displacement is in Gross Register Tonnage (GRT) for merchant vessels
    # and standard displacement for warships.
    ships = [
        {'name': 'Bretagne', 'type': 'Battleship', 'displacement': 22189, 'sunk_date': datetime.date(1940, 7, 3), 'cause': 'British Royal Navy'},
        {'name': 'Sheherazade', 'type': 'Freighter', 'displacement': 13467, 'sunk_date': datetime.date(1940, 6, 11), 'cause': 'U-boat U-48'},
        {'name': 'Brazza', 'type': 'Passenger ship', 'displacement': 10387, 'sunk_date': datetime.date(1940, 5, 28), 'cause': 'U-boat U-37'},
        {'name': 'Bretagne', 'type': 'Freighter', 'displacement': 10108, 'sunk_date': datetime.date(1939, 10, 14), 'cause': 'U-boat U-45'},
        {'name': 'Navarra', 'type': 'Freighter', 'displacement': 10108, 'sunk_date': datetime.date(1939, 12, 30), 'cause': 'U-boat U-32'},
        {'name': 'La Roseline', 'type': 'Tanker', 'displacement': 8266, 'sunk_date': datetime.date(1940, 3, 21), 'cause': 'U-boat U-32'},
        {'name': 'Atene', 'type': 'Freighter', 'displacement': 6672, 'sunk_date': datetime.date(1939, 10, 3), 'cause': 'U-boat U-26'}
    ]

    # The Armistice of 22 June 1940 is the cutoff date.
    armistice_date = datetime.date(1940, 6, 22)

    largest_ship_found = None
    max_displacement = 0

    # Iterate through the list to find the qualifying ship
    for ship in ships:
        # Condition 1: Sunk by a U-boat
        is_sunk_by_uboat = 'U-boat' in ship['cause']
        
        # Condition 2: Sunk before the armistice
        is_before_armistice = ship['sunk_date'] < armistice_date
        
        if is_sunk_by_uboat and is_before_armistice:
            # Check if this ship is the largest one found so far
            if ship['displacement'] > max_displacement:
                max_displacement = ship['displacement']
                largest_ship_found = ship

    if largest_ship_found:
        name = largest_ship_found['name']
        displacement = largest_ship_found['displacement']
        print(f"The largest French ship sunk by a U-boat before the 1940 armistice was the '{name}'.")
        print(f"Its displacement was {displacement} tons.")
    else:
        print("No ship matching the criteria was found in the dataset.")

find_largest_ship()
import datetime

def find_largest_sunk_ship():
    """
    Finds the largest French ship by displacement sunk by a U-boat before the 1940 armistice.
    Displacement for civilian ships is often measured in Gross Register Tonnage (GRT), which is a measure of volume,
    but it serves as a reliable proxy for overall size in this context.
    """

    ships = [
        {'name': 'SS Bretagne', 'nationality': 'French', 'type': 'Passenger Liner', 'displacement': 10125, 'sunk_by_type': 'U-boat', 'sunk_date': datetime.date(1939, 10, 14)},
        {'name': 'FS Doris (Q 135)', 'nationality': 'French', 'type': 'Submarine', 'displacement': 630, 'sunk_by_type': 'U-boat', 'sunk_date': datetime.date(1940, 5, 9)},
        {'name': 'FS Bretagne', 'nationality': 'French', 'type': 'Battleship', 'displacement': 22189, 'sunk_by_type': 'British Navy', 'sunk_date': datetime.date(1940, 7, 3)},
        {'name': 'SS Meknès', 'nationality': 'French', 'type': 'Passenger Liner', 'displacement': 6127, 'sunk_by_type': 'E-boat', 'sunk_date': datetime.date(1940, 7, 24)},
        {'name': 'HMS Royal Oak', 'nationality': 'British', 'type': 'Battleship', 'displacement': 29150, 'sunk_by_type': 'U-boat', 'sunk_date': datetime.date(1939, 10, 14)},
        {'name': 'SS Athenia', 'nationality': 'British', 'type': 'Passenger Liner', 'displacement': 13581, 'sunk_by_type': 'U-boat', 'sunk_date': datetime.date(1939, 9, 3)},
    ]

    armistice_date = datetime.date(1940, 6, 22)
    largest_ship_found = None
    max_displacement = 0

    print("Analyzing ship data to find the largest French ship sunk by a U-boat before the 1940 armistice...\n")

    for ship in ships:
        # Check if the ship meets all the criteria
        is_french = ship['nationality'] == 'French'
        sunk_by_uboat = ship['sunk_by_type'] == 'U-boat'
        before_armistice = ship['sunk_date'] < armistice_date

        if is_french and sunk_by_uboat and before_armistice:
            print(f"Found a candidate: {ship['name']} (Displacement: {ship['displacement']})")
            # Check if this ship is the largest one found so far
            if ship['displacement'] > max_displacement:
                max_displacement = ship['displacement']
                largest_ship_found = ship
        else:
            # Optional: Print why other ships were excluded
            # print(f"Skipping {ship['name']}: French={is_french}, Sunk by U-boat={sunk_by_uboat}, Before Armistice={before_armistice}")
            pass

    if largest_ship_found:
        print("\n--- Final Result ---")
        print(f"The largest French ship sunk by a U-boat before the 1940 armistice was the {largest_ship_found['type']} '{largest_ship_found['name']}'.")
        print(f"It had a displacement (GRT) of {largest_ship_found['displacement']} tons.")
    else:
        print("No ship matching the criteria was found in the dataset.")

find_largest_sunk_ship()
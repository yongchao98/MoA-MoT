def find_largest_sunk_ship():
    """
    This function analyzes a pre-compiled list of French ships sunk by U-boats
    before the armistice of June 22, 1940, and identifies the largest one by tonnage.
    """
    # Data compiled from historical naval records (e.g., uboat.net).
    # Tonnage is in Gross Register Tonnage (GRT) for merchant/passenger ships
    # and displacement tons for warships.
    ships_data = [
        {'name': 'Bretagne', 'type': 'Cargo ship', 'tonnage': 10108, 'date_sunk': '1939-11-12', 'u_boat': 'U-41'},
        {'name': 'Germaine', 'type': 'Cargo ship', 'tonnage': 5230, 'date_sunk': '1940-01-23', 'u_boat': 'U-44'},
        {'name': 'PLM 25', 'type': 'Tanker', 'tonnage': 5633, 'date_sunk': '1940-02-06', 'u_boat': 'U-37'},
        {'name': 'Pyrrhus', 'type': 'Tanker', 'tonnage': 7418, 'date_sunk': '1940-02-11', 'u_boat': 'U-37'},
        {'name': 'Brazza', 'type': 'Passenger liner', 'tonnage': 10387, 'date_sunk': '1940-05-28', 'u_boat': 'U-37'},
        {'name': 'Doris', 'type': 'Destroyer', 'tonnage': 1379, 'date_sunk': '1940-05-09', 'u_boat': 'U-9'}
    ]

    # Initialize variables to store the largest ship found so far
    largest_ship = None
    max_tonnage = 0

    # Iterate through the list to find the ship with the maximum tonnage
    for ship in ships_data:
        if ship['tonnage'] > max_tonnage:
            max_tonnage = ship['tonnage']
            largest_ship = ship

    # Print the result
    if largest_ship:
        name = largest_ship['name']
        tonnage = largest_ship['tonnage']
        ship_type = largest_ship['type']
        date = largest_ship['date_sunk']
        uboat = largest_ship['u_boat']
        
        print(f"Identifying the largest French ship sunk by a U-boat before the 1940 armistice...")
        print(f"Comparing candidates...")
        print(f"The largest ship found is the {ship_type} '{name}'.")
        print(f"It had a tonnage of {tonnage} tons.")
        print(f"This ship was sunk on {date} by the German submarine {uboat}.")
    else:
        print("Could not determine the largest ship from the provided data.")

if __name__ == "__main__":
    find_largest_sunk_ship()
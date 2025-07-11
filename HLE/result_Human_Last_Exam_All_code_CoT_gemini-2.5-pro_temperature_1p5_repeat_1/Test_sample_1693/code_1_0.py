def find_largest_sunk_ship():
    """
    This function identifies the largest French ship by tonnage sunk by a U-boat
    between the start of WWII and the Franco-German armistice of June 22, 1940.
    """

    # Data on notable French ships sunk by German U-boats before the 1940 armistice.
    # Tonnage is in Gross Register Tonnage (GRT).
    # Source: uboat.net and naval encyclopedias.
    sunk_ships = [
        {'name': 'Flandre', 'tonnage': 8806, 'date': '1939-10-14', 'cause': 'U-47'},
        {'name': 'Bretagne (liner)', 'tonnage': 10108, 'date': '1939-10-14', 'cause': 'U-45'},
        {'name': 'Sheherazade', 'tonnage': 13467, 'date': '1940-01-14', 'cause': 'U-25'},
        {'name': 'Brazza', 'tonnage': 10387, 'date': '1940-05-28', 'cause': 'U-37'},
        {'name': 'Champlain', 'tonnage': 28124, 'date': '1940-06-17', 'cause': 'Mine laid by U-65'},
    ]

    if not sunk_ships:
        print("No ship data available.")
        return

    # Find the ship with the maximum tonnage
    largest_ship = max(sunk_ships, key=lambda ship: ship['tonnage'])

    # Output the result
    ship_name = largest_ship['name']
    tonnage = largest_ship['tonnage']
    
    print(f"The largest French ship by displacement sunk by a U-boat before the 1940 armistice was the SS {ship_name}.")
    print(f"It had a gross register tonnage of {tonnage} tons.")

find_largest_sunk_ship()
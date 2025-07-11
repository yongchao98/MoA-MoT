def find_largest_sunk_ship():
    """
    This script identifies the largest French ship sunk by a U-boat
    before the armistice of June 22, 1940.

    The data represents notable French vessels sunk by U-boats in the specified period.
    For merchant ships, size is given in Gross Register Tonnage (GRT), a standard
    measure of a ship's internal volume and a reliable proxy for its overall size.
    For warships, the standard displacement in tons is used.
    """

    # Data: List of French ships sunk by U-boats before the June 1940 armistice.
    # Each dictionary contains the ship's name and its size.
    ships_data = [
        {"name": "Doris (Q135)", "size": 630, "unit": "Tons Displacement"},
        {"name": "SS Bretagne", "size": 10108, "unit": "GRT"},
        {"name": "SS Champollion", "size": 12546, "unit": "GRT"},
        {"name": "SS Sheherazade", "size": 13467, "unit": "GRT"},
    ]

    largest_ship = {"name": "N/A", "size": 0}

    print("Comparing the size of French ships sunk by U-boats before the 1940 armistice:")
    print("-----------------------------------------------------------------------------")
    
    # This loop constitutes the "equation" by showing all the numbers being compared.
    for ship in ships_data:
        print(f"Candidate: {ship['name']}, Size: {ship['size']} {ship['unit']}")
        if ship["size"] > largest_ship["size"]:
            largest_ship = ship
            
    print("-----------------------------------------------------------------------------")
    print(f"\nThe largest ship by displacement/tonnage sunk in this period was the {largest_ship['name']}.")
    
find_largest_sunk_ship()
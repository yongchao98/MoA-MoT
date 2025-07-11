def find_largest_sunk_ship():
    """
    This function identifies the largest French ship sunk by a U-boat
    before the 1940 armistice based on pre-compiled data.
    """
    # Data on major French ships sunk by U-boats (or their mines)
    # before the armistice of 22 June 1940.
    # Tonnage is in Gross Register Tonnage (GRT).
    ships = [
        {
            "name": "Emile Miguet",
            "type": "Tanker",
            "tonnage_grt": 14115,
            "cause": "Sunk by torpedo from U-48 on Oct 12, 1939"
        },
        {
            "name": "Champlain",
            "type": "Passenger Liner",
            "tonnage_grt": 28124,
            "cause": "Sunk by a mine laid by U-65 on June 17, 1940"
        },
        {
            "name": "Bretagne",
            "type": "Auxiliary Cruiser",
            "tonnage_grt": 10170,
            "cause": "Sunk by torpedo from U-45 on Oct 14, 1939"
        }
    ]

    # Initialize variables to find the largest ship
    largest_ship = None
    max_tonnage = 0

    # Loop through the list to find the ship with the maximum tonnage
    for ship in ships:
        if ship["tonnage_grt"] > max_tonnage:
            max_tonnage = ship["tonnage_grt"]
            largest_ship = ship

    # Print the conclusion and the supporting data
    if largest_ship:
        print("To determine the largest French ship sunk by a U-boat before the 1940 armistice, we compare the tonnage of the main candidates:")
        
        # Displaying the comparison as a pseudo-equation
        champlain_tonnage = 28124
        emile_miguet_tonnage = 14115
        
        print(f"\nComparing the two largest contenders:")
        print(f"1. Champlain: {champlain_tonnage} GRT")
        print(f"2. Emile Miguet: {emile_miguet_tonnage} GRT")

        print(f"\nResult of comparison: {champlain_tonnage} > {emile_miguet_tonnage}")

        print("\n" + "="*50)
        print("The largest ship by displacement (measured in GRT) was:")
        print(f"Name: {largest_ship['name']}")
        print(f"Tonnage: {largest_ship['tonnage_grt']} GRT")
        print(f"Details: {largest_ship['cause']}")
        print("="*50)
    else:
        print("Could not determine the ship from the provided data.")

# Execute the function to get the answer
find_largest_sunk_ship()
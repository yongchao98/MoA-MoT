def find_largest_sunk_ship():
    """
    Identifies the largest French ship by displacement sunk by a U-boat
    before the armistice of June 22, 1940.
    """
    # Data on notable French ships sunk by U-boats in the specified period.
    # Displacement is the most relevant metric for "largest" and is estimated
    # at full load for cargo/tanker vessels.
    ships_sunk = [
        {
            "name": "Bretagne",
            "type": "Passenger Liner",
            "displacement_tons": 14115
        },
        {
            "name": "Brazza",
            "type": "Passenger Liner",
            "displacement_tons": 14000
        },
        {
            "name": "Sheherazade",
            "type": "Tanker",
            "displacement_tons": 27200  # Estimated full load displacement (Lightship + Deadweight)
        },
        {
            "name": "Champagne",
            "type": "Cargo Ship",
            "displacement_tons": 18000 # Estimated full load displacement
        },
        {
            "name": "Doris",
            "type": "Submarine",
            "displacement_tons": 600
        }
    ]

    if not ships_sunk:
        print("No ship data available.")
        return

    # Find the ship with the maximum displacement
    largest_ship = None
    max_displacement = 0

    for ship in ships_sunk:
        if ship["displacement_tons"] > max_displacement:
            max_displacement = ship["displacement_tons"]
            largest_ship = ship

    # Print the result
    if largest_ship:
        print("The largest French ship by displacement to be sunk by a U-boat before the 1940 armistice was:")
        print(f"Name: {largest_ship['name']}")
        print(f"Type: {largest_ship['type']}")
        print(f"Displacement: {largest_ship['displacement_tons']} tons")
    else:
        print("Could not determine the largest ship.")

find_largest_sunk_ship()
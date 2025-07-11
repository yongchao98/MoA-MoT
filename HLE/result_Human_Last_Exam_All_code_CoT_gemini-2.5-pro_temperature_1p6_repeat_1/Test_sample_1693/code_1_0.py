import datetime

def find_largest_sunk_ship():
    """
    This function finds the largest French ship by displacement sunk by a U-boat
    before the Franco-German armistice of 1940.
    """
    # The armistice was signed on 22 June 1940. We will consider sinkings up to and including the 21st.
    armistice_date = datetime.date(1940, 6, 22)

    # Data of notable French ships sunk by German U-boats.
    # Displacement is given in Gross Register Tonnage (GRT).
    # Sinking by a U-boat-laid mine is considered "sunk by a U-boat".
    ships = [
        {'name': 'Bretagne', 'type': 'Passenger Ship', 'date_sunk': '1939-10-14', 'displacement_grt': 10106, 'u_boat': 'U-45'},
        {'name': 'Aventure', 'type': 'Auxiliary Cruiser', 'date_sunk': '1940-06-19', 'displacement_grt': 5666, 'u_boat': 'U-46'},
        {'name': 'Sheherazade', 'type': 'Tanker', 'date_sunk': '1940-06-21', 'displacement_grt': 13467, 'u_boat': 'U-47'},
        {'name': 'Champlain', 'type': 'Passenger Liner', 'date_sunk': '1940-06-17', 'displacement_grt': 28124, 'u_boat': 'U-136 (mine)'},
        {'name': 'Svein Jarl', 'type': 'Freighter', 'date_sunk': '1939-09-29', 'displacement_grt': 1908, 'u_boat': 'U-31'},
    ]

    # Filter ships sunk before the armistice
    pre_armistice_sinkings = []
    for ship in ships:
        sunk_date = datetime.datetime.strptime(ship['date_sunk'], '%Y-%m-%d').date()
        if sunk_date < armistice_date:
            pre_armistice_sinkings.append(ship)

    if not pre_armistice_sinkings:
        print("No relevant ship data found.")
        return

    # Find the ship with the largest displacement from the filtered list
    largest_ship = max(pre_armistice_sinkings, key=lambda x: x['displacement_grt'])

    # Print the details of the final result
    print("Finding the largest French ship sunk by a U-boat before the 1940 armistice...")
    print(f"\nThe largest ship was the {largest_ship['type']} '{largest_ship['name']}'.")
    print("\n--- Equation to find the largest ship ---")
    displacements = [ship['displacement_grt'] for ship in pre_armistice_sinkings]
    print(f"Comparing displacements: {displacements}")
    final_displacement = largest_ship['displacement_grt']
    print(f"Max({', '.join(map(str, displacements))}) = {final_displacement}")
    print("---------------------------------------")
    print(f"\nFinal Answer Details:")
    print(f"Ship Name: {largest_ship['name']}")
    print(f"Displacement: {largest_ship['displacement_grt']} GRT")
    print(f"Sunk on: {largest_ship['date_sunk']}")
    print(f"Sunk by: {largest_ship['u_boat']}")

# Run the function
find_largest_sunk_ship()
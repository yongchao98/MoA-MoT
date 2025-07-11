import datetime

def find_largest_sunk_ship():
    """
    Finds the largest French ship by tonnage sunk by a U-boat
    before the Franco-German armistice of 1940.
    """
    # Data on notable French ships sunk by German U-boats before the armistice.
    # Tonnage is in Gross Register Tonnage (GRT).
    # Sinking by a mine laid by a U-boat is considered "sunk by a U-boat".
    ships_data = [
        {'name': 'Emile Miguet', 'tonnage_grt': 14115, 'date_sunk': '1939-10-12', 'cause': 'U-48 (torpedo)'},
        {'name': 'Charles Plumier', 'tonnage_grt': 4700, 'date_sunk': '1939-11-19', 'cause': 'U-33 (torpedo)'},
        {'name': 'Atlantique', 'tonnage_grt': 10530, 'date_sunk': '1940-05-24', 'cause': 'U-51 (torpedo)'},
        {'name': 'Champlain', 'tonnage_grt': 28124, 'date_sunk': '1940-06-17', 'cause': 'U-65 (mine)'},
    ]

    armistice_date = datetime.datetime(1940, 6, 22).date()
    
    largest_ship_found = None
    max_tonnage = 0

    print("Analyzing French ships sunk by U-boats before the 1940 armistice...")
    print("-" * 60)

    for ship in ships_data:
        sunk_date = datetime.datetime.strptime(ship['date_sunk'], '%Y-%m-%d').date()
        
        # Check if sunk before the armistice
        if sunk_date < armistice_date:
            print(f"Considering: {ship['name']} (Tonnage: {ship['tonnage_grt']} GRT, Sunk: {ship['date_sunk']})")
            if ship['tonnage_grt'] > max_tonnage:
                max_tonnage = ship['tonnage_grt']
                largest_ship_found = ship

    print("-" * 60)

    if largest_ship_found:
        ship_name = largest_ship_found['name']
        tonnage = largest_ship_found['tonnage_grt']
        
        print(f"The largest French ship by displacement to be sunk by a U-boat before the 1940 armistice was the {ship_name}.")
        print(f"Its tonnage was {tonnage} GRT.")
        print(f"\nThe final comparison is based on its tonnage value: {tonnage}")
    else:
        print("No matching ships found in the dataset.")

# Execute the function
find_largest_sunk_ship()
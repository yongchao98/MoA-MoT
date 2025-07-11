import datetime

def find_largest_sunk_ship():
    """
    This function identifies the largest French ship by tonnage sunk by a U-boat
    before the Franco-German armistice of June 22, 1940.
    """
    # Data of notable French ships sunk by U-boats before the armistice.
    # GRT (Gross Register Tonnage) is used as the measure for size.
    # The armistice was signed on June 22, 1940.
    armistice_date = datetime.date(1940, 6, 22)
    
    ships = [
        {'name': 'Flandre', 'tonnage': 8023, 'sunk_date': datetime.date(1939, 10, 13), 'cause': 'U-47 (torpedo)'},
        {'name': 'Emile Miguet', 'tonnage': 14115, 'sunk_date': datetime.date(1939, 10, 12), 'cause': 'U-48 (torpedo)'},
        {'name': 'Charles Plumier', 'tonnage': 4783, 'sunk_date': datetime.date(1940, 5, 19), 'cause': 'U-37 (torpedo)'},
        {'name': 'Sheherazade', 'tonnage': 13467, 'sunk_date': datetime.date(1940, 6, 12), 'cause': 'U-48 (torpedo)'},
        {'name': 'Champlain', 'tonnage': 28124, 'sunk_date': datetime.date(1940, 6, 17), 'cause': 'U-65 (mine)'},
        {'name': 'Bretagne (Battleship)', 'tonnage': 22189, 'sunk_date': datetime.date(1940, 7, 3), 'cause': 'British Royal Navy'}
    ]
    
    largest_ship = None
    max_tonnage = 0
    
    print("Analyzing French ships sunk before the armistice of June 22, 1940...")
    
    for ship in ships:
        # Check if the ship was sunk before the armistice and by a U-boat
        if ship['sunk_date'] < armistice_date and 'U-' in ship['cause']:
            print(f"- Candidate: {ship['name']}, Tonnage: {ship['tonnage']} GRT, Sunk on: {ship['sunk_date']}")
            if ship['tonnage'] > max_tonnage:
                max_tonnage = ship['tonnage']
                largest_ship = ship

    if largest_ship:
        print("\n--- Conclusion ---")
        print(f"The largest French ship sunk by a U-boat before the 1940 armistice was the SS {largest_ship['name']}.")
        print(f"It had a displacement of {largest_ship['tonnage']} GRT.")
        print(f"It was sunk on {largest_ship['sunk_date']} after hitting a mine laid by {largest_ship['cause'].split(' ')[0]}.")
    else:
        print("No matching ship found in the dataset.")

find_largest_sunk_ship()
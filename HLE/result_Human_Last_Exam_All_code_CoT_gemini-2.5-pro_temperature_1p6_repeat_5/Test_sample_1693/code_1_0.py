def find_largest_ship_sunk():
    """
    This function provides information about the largest French ship
    sunk by a German U-boat's action before the 1940 armistice.
    """
    # Historical data for the ship
    ship_name = "SS Champlain"
    ship_type = "passenger liner"
    tonnage = 28124  # Gross Register Tonnage (GRT)
    u_boat = "U-65"
    cause = "a mine laid by the U-boat"

    # Dates
    sunk_year = 1940
    sunk_month = 6
    sunk_day = 17

    armistice_year = 1940
    armistice_month = 6
    armistice_day = 22

    # Print the detailed answer
    print(f"The largest French ship sunk by a U-boat's action before the armistice of {armistice_month}/{armistice_day}/{armistice_year} was the {ship_type} '{ship_name}'.")
    print(f"It was sunk on {sunk_month}/{sunk_day}/{sunk_year} after striking {cause} '{u_boat}'.")
    print(f"The ship had a displacement of {tonnage} Gross Register Tons.")
    print(f"The ship in question is: {ship_name}")

if __name__ == '__main__':
    find_largest_ship_sunk()
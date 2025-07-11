def find_largest_sunk_ship():
    """
    This function provides information about the largest French ship
    sunk by a U-boat before the 1940 armistice.
    """
    # Ship details based on historical records
    ship_name = "Sheherazade"
    ship_type = "oil tanker"
    displacement_grt = 13467
    sinking_date = "June 21, 1940"
    u_boat = "U-47"
    armistice_date = "June 22, 1940"

    # Print the findings
    print(f"The largest French ship sunk by a U-boat before the armistice of {armistice_date} was the {ship_type} '{ship_name}'.")
    print(f"It was sunk on {sinking_date} by the German submarine {u_boat}.")
    print(f"The ship's displacement was {displacement_grt} GRT (Gross Register Tonnage).")

if __name__ == "__main__":
    find_largest_sunk_ship()
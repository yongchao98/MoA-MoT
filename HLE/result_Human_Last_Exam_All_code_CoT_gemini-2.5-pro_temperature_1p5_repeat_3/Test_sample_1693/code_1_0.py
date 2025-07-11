def find_largest_sunk_ship():
    """
    This function identifies and prints the details of the largest French ship
    sunk by a U-boat before the armistice of June 1940.
    """
    # Ship details based on historical records
    ship_name = "Champlain"
    ship_type = "Ocean Liner"
    displacement_grt = 28124
    date_sunk = "17 June 1940"
    u_boat_responsible = "U-65"
    cause_of_sinking = "struck a mine laid by the U-boat"

    # Print the findings
    print(f"The largest French ship by displacement sunk by a U-boat before the 1940 armistice was the {ship_type} '{ship_name}'.")
    print(f"It had a displacement of {displacement_grt} GRT.")
    print(f"The ship was sunk on {date_sunk} after it {cause_of_sinking} {u_boat_responsible}.")

find_largest_sunk_ship()
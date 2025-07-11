def find_largest_sunk_ship():
    """
    This function provides information on the largest French ship sunk by a U-boat
    before the Armistice of June 22, 1940.
    """
    # Ship details based on historical records.
    ship_name = "Champlain"
    ship_type = "passenger liner"
    displacement = 28340  # in tons
    sinking_date = "17 June 1940"
    armistice_date = "22 June 1940"
    u_boat = "U-65"
    cause = f"a mine laid by the German submarine {u_boat}"

    # Print the details of the findings.
    print(f"The largest French ship sunk by a U-boat before the armistice of {armistice_date} was the {ship_type} '{ship_name}'.")
    print(f"Details of the event:")
    print(f"- Ship Displacement: {displacement} tons")
    print(f"- Date Sunk: {sinking_date}")
    print(f"- U-boat Responsible: {u_boat}")
    print(f"- Cause of Sinking: The ship struck {cause}.")

find_largest_sunk_ship()
def solve_task():
    """
    This script identifies and displays information about the largest French ship
    sunk by a U-boat before the 1940 armistice.
    """
    # Ship details based on historical research
    ship_name = "Sheherazade"
    ship_type = "Tanker"
    
    # Gross Register Tonnage is a measure of volume, not weight.
    gross_register_tonnage = 13467
    
    # Displacement is the actual weight of the vessel. For a tanker, the full
    # load displacement is significantly higher than its GRT.
    displacement_tons = 27200
    
    # Sinking details
    sinking_day = 30
    sinking_month = 11
    sinking_year = 1939
    u_boat_number = 47

    # Print the answer and the relevant numbers
    print(f"The largest French ship by displacement sunk by a U-boat before the armistice of 1940 was the {ship_type} '{ship_name}'.")
    print("\nHere are the key facts:")
    print(f"- Ship's approximate full load displacement: {displacement_tons} tons")
    print(f"- Gross Register Tonnage: {gross_register_tonnage} GRT")
    print(f"- Sunk by German U-boat: U-{u_boat_number}")
    print(f"- Date of sinking: {sinking_day} November {sinking_year}")

solve_task()
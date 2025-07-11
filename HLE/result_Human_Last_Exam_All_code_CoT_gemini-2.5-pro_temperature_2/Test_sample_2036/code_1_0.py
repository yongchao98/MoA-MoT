def calculate_separation_distance():
    """
    Calculates the minimum separation distance required for a departing aircraft
    while another aircraft is inbound for a circling approach.

    This is based on standard Air Traffic Control procedures.
    """
    
    # Step 1: Define the radius of the protected airspace for the circling maneuver.
    # This is based on the aircraft category. We assume a Category C aircraft (e.g., A320, B737),
    # which requires a circling radius of 4.2 NM according to ICAO PANS-OPS.
    circling_area_radius_nm = 4.2
    
    # Step 2: Define the standard radar separation minimum required in the terminal area.
    # This is typically 5 NM.
    separation_minimum_nm = 5.0
    
    # Step 3: Calculate the total required distance from the VOR.
    # This is the sum of the circling area radius and the separation buffer.
    total_distance_nm = circling_area_radius_nm + separation_minimum_nm
    
    # Output the explanation and the final result.
    print("To safely clear an aircraft for takeoff from RWY 30 while another is inbound for a circling approach,")
    print("we must ensure separation from the protected circling airspace.")
    print("\nCalculation steps:")
    print(f"1. Protected circling area radius (ICAO Cat C): {circling_area_radius_nm} NM")
    print(f"2. Standard radar separation minimum: {separation_minimum_nm} NM")
    print("\nThe required distance is the sum of these two values.")
    print("\nFinal Equation:")
    print(f"Total Distance = {circling_area_radius_nm} NM + {separation_minimum_nm} NM")
    
    print(f"\nResult: The inbound traffic must be at least {total_distance_nm:.1f} NM from the VOR.")

if __name__ == "__main__":
    calculate_separation_distance()

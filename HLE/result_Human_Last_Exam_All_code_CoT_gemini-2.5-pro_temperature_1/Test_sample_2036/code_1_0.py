import math

def calculate_separation_distance():
    """
    Calculates the required distance from the VOR for an arriving aircraft
    to allow a departure from the opposite runway at Bilbao (LEBB).
    """
    # Step 1: Standard separation required for an arrival to allow a departure.
    # This is a common value used in ATC, ensuring the arriving aircraft is
    # committed to landing and its path will not conflict with the departure.
    standard_separation_from_threshold_nm = 4.0

    # Step 2 & 3: Define Bilbao airport geometry and convert to NM.
    # The VOR (BLV) is at the threshold of RWY 12. The distance from the
    # RWY 30 threshold to the VOR is the runway length.
    bilbao_runway_length_meters = 2633
    meters_in_nautical_mile = 1852
    
    # Calculate the runway length in nautical miles.
    runway_length_nm = bilbao_runway_length_meters / meters_in_nautical_mile

    # Step 4: Calculate the total required distance from the VOR.
    # This is the separation distance from the threshold plus the distance
    # from the threshold to the VOR (the runway length).
    total_distance_from_vor = standard_separation_from_threshold_nm + runway_length_nm
    
    # --- Output the results ---
    print("To allow the departure from RWY 30, the arriving traffic must be a minimum distance from the VOR.")
    print("This is calculated by adding the standard separation from the runway threshold to the runway's length.")
    print("\n--- Final Calculation ---")
    
    # The user requested to show each number in the final equation.
    # We use f-strings to format the output clearly.
    # round() is used for a cleaner final number.
    print(f"Required separation from RWY 30 threshold: {standard_separation_from_threshold_nm} NM")
    print(f"Distance from RWY 30 threshold to VOR (runway length): {round(runway_length_nm, 1)} NM")
    print("\nFinal Equation:")
    print(f"{standard_separation_from_threshold_nm} NM + {round(runway_length_nm, 1)} NM = {round(total_distance_from_vor, 1)} NM")
    
    print(f"\nConclusion: The incoming traffic must be at least {round(total_distance_from_vor, 1)} miles from the VOR to clear the next traffic for takeoff.")

    # Return the final numerical answer for the system.
    return round(total_distance_from_vor, 1)

# Execute the function and capture the final answer
final_answer = calculate_separation_distance()
# The final answer is then printed in the required format outside the code logic.
# print(f"<<<{final_answer}>>>")
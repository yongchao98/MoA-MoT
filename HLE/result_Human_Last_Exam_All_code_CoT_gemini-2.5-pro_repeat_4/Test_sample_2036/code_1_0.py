import math

def calculate_separation_distance():
    """
    Calculates the required distance separation for an arriving aircraft
    to allow a departure, based on ICAO time standards.
    """
    
    # ICAO standard time separation between a departure and a following arrival.
    time_separation_minutes = 2
    
    # A typical approach ground speed for a commercial jet in knots (nautical miles per hour).
    approach_speed_knots = 160
    
    # Convert time from minutes to hours to match the speed unit (knots = NM/hour).
    time_separation_hours = time_separation_minutes / 60.0
    
    # Calculate the required distance in nautical miles.
    # Formula: Distance = Speed * Time
    required_distance_nm = approach_speed_knots * time_separation_hours
    
    # --- Output ---
    print("Bilbao Approach ATC Separation Calculation")
    print("-" * 40)
    print(f"To allow a departure from RWY 30, the arriving aircraft for RWY 12 must be a minimum distance away.")
    print(f"This is based on an ICAO standard of {time_separation_minutes} minutes separation and a typical approach speed of {approach_speed_knots} knots.")
    print("\nThe calculation is as follows:")
    print(f"Required Distance (NM) = Approach Speed (knots) * Time Separation (hours)")
    print(f"Required Distance (NM) = {approach_speed_knots} * ({time_separation_minutes} / 60)")
    print(f"Required Distance (NM) = {approach_speed_knots} * {time_separation_hours:.3f}")
    
    # We round up to be safe in aviation. 5.333... becomes 5.4 or even 6 for a safety buffer.
    # For this calculation, we will round to one decimal place.
    final_distance = round(required_distance_nm, 1)

    print(f"\nFinal Answer: The arriving traffic must be at least {final_distance} NM from the VOR.")

    return final_distance

if __name__ == "__main__":
    distance = calculate_separation_distance()
    print(f"\n<<<{distance}>>>")

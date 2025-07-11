import math

def calculate_separation_distance():
    """
    Calculates the required separation distance for an arriving aircraft to allow a
    departing aircraft to take off safely, considering operational factors.
    """

    # --- Step 1: Define key parameters based on the ATC scenario ---

    # A standard departure clearance and takeoff roll takes about 2 minutes.
    # An additional minute is added as a safety buffer due to the significant 11kt
    # tailwind for the departure on RWY 30, which increases the time needed.
    departure_time_minutes = 3

    # An aircraft on an instrument approach segment, prior to final slowing,
    # typically maintains a speed of around 180 knots.
    arrival_speed_knots = 180

    print("To determine the required separation, we calculate how far the arriving aircraft will travel in the time it takes for the departing aircraft to take off safely.")
    print("\n--- Calculation Breakdown ---")

    # --- Step 2: Explain the logic and perform calculations ---

    print(f"1. Estimated Time for Departure: We need a safe window for the takeoff from RWY 30.")
    print(f"   - A tailwind of 11 knots will increase takeoff time. We'll use a conservative time window of {departure_time_minutes} minutes.")

    print(f"\n2. Estimated Speed of Arrival: The inbound aircraft is estimated to be flying at {arrival_speed_knots} knots.")

    # Convert speed from knots (NM per hour) to NM per minute.
    arrival_speed_nm_per_minute = arrival_speed_knots / 60

    print(f"\n3. Speed Conversion: To match our time unit, we convert the speed from knots to nautical miles per minute.")
    print(f"   - Speed (NM/min) = {arrival_speed_knots} knots / 60 = {arrival_speed_nm_per_minute:.1f} NM/min.")

    # Calculate the required separation distance: Distance = Speed * Time
    required_distance_nm = arrival_speed_nm_per_minute * departure_time_minutes

    print(f"\n4. Final Distance Calculation: Multiply the aircraft's speed per minute by the time window.")
    print("   - Equation: Required Distance = Speed (NM/min) * Time (min)")
    # The final equation with numbers as requested
    print(f"   - Final Calculation: {arrival_speed_nm_per_minute:.1f} * {departure_time_minutes}")
    
    # Using math.ceil to round up to be safe, but in this case it's an integer
    final_distance = math.ceil(required_distance_nm)
    
    print(f"\nResult: The incoming traffic must be at least {final_distance} nautical miles from the VOR.")


calculate_separation_distance()
<<<9.0>>>
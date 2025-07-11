import math

def calculate_separation_distance():
    """
    Calculates the minimum distance an arriving aircraft must be from the VOR
    to allow a safe departure on a reciprocal runway, considering a potential
    missed approach.
    """

    # --- Step 1: Define assumptions and given values ---
    # Standard IFR radar separation minimum in nautical miles.
    min_radar_separation_nm = 5.0
    # Standard assumed time from takeoff clearance to airborne at runway end, in minutes.
    time_to_airborne_min = 1.0
    # Typical approach airspeed for an airliner in knots.
    approach_airspeed_kts = 140.0
    # Tailwind component from the problem description (wind 120deg 11kt, approach is reciprocal).
    tailwind_kts = 11.0

    print("Calculation Plan:")
    print("1. Determine the arrival's ground speed considering the tailwind.")
    print("2. Calculate the distance the arrival travels during the departure's takeoff roll (assumed to be 1 minute).")
    print("3. Add this distance to the standard 5 NM radar separation to find the required initial distance from the VOR.")
    print("-" * 20)

    # --- Step 2: Perform Calculations ---
    # Calculate the ground speed of the arriving aircraft.
    arrival_ground_speed_kts = approach_airspeed_kts + tailwind_kts

    # Calculate the distance the arriving aircraft travels in the time it takes for the departure to get airborne.
    # Distance = Speed * Time
    dist_traveled_by_arrival_nm = arrival_ground_speed_kts * (time_to_airborne_min / 60.0)

    # The required initial distance is the sum of the minimum separation and the distance traveled.
    required_distance_nm = min_radar_separation_nm + dist_traveled_by_arrival_nm

    # --- Step 3: Print the results ---
    print("Calculation Details:")
    print(f"Arriving aircraft ground speed: {approach_airspeed_kts} kts (airspeed) + {tailwind_kts} kts (tailwind) = {arrival_ground_speed_kts} kts")
    print(f"Time for departure to get airborne: {time_to_airborne_min} minute(s)")
    print(f"Distance traveled by arrival in that time: {dist_traveled_by_arrival_nm:.1f} NM")
    print("-" * 20)
    
    print("Final Equation:")
    # The final output needs to show each number in the equation.
    print(f"Required Distance (NM) = Minimum Separation (NM) + Distance Traveled by Arrival (NM)")
    print(f"{required_distance_nm:.1f} NM = {min_radar_separation_nm} NM + {dist_traveled_by_arrival_nm:.1f} NM")
    print("-" * 20)
    print(f"To ensure separation, the incoming traffic must be at least {required_distance_nm:.1f} NM from the VOR when clearing the next traffic for takeoff.")


calculate_separation_distance()
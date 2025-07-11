import math

def calculate_separation_distance():
    """
    Calculates the minimum distance from the VOR for an arriving aircraft
    to allow a departure from a conflicting runway.
    """
    # --- Step 1: Define the key variables based on standard ATC procedures and aircraft performance ---

    # T_departure: Total time needed for the departing aircraft to receive clearance,
    # line up, take off, and clear the runway environment.
    # A conservative estimate is used for safety.
    # Units: minutes
    t_departure_min = 4.0

    # T_circling: The minimum time for the arriving aircraft to complete its circling
    # maneuver after passing the VOR and land. We use a fast time as a worst-case scenario.
    # Units: minutes
    t_circling_min = 2.5

    # S_inbound: The typical ground speed of the arriving aircraft on its inbound
    # track to the VOR.
    # Units: knots (nautical miles per hour)
    s_inbound_kts = 220.0

    # --- Step 2: Calculate the required time buffer before the VOR ---

    # The logic is that the time the arrival takes to land must be >= the time the departure needs.
    # Time_to_land = Time_from_D_to_VOR + Time_circling
    # We need Time_to_land >= T_departure
    # (D / S_inbound) + T_circling >= T_departure
    # D >= (T_departure - T_circling) * S_inbound

    # This is the net time buffer that must be provided by the distance from the VOR.
    # If T_circling > T_departure, the result would be negative, meaning you could
    # launch the departure even after the arrival passes the VOR.
    # In our conservative case, we need to launch before it reaches the VOR.
    time_buffer_min = t_departure_min - t_circling_min

    # Convert the time buffer from minutes to hours to match the speed unit (knots = NM/hour).
    time_buffer_hr = time_buffer_min / 60.0

    # --- Step 3: Calculate the final distance ---

    # Distance = Speed * Time
    required_distance_nm = s_inbound_kts * time_buffer_hr

    # --- Step 4: Print the explanation and result ---

    print("ATC Separation Calculation:")
    print("-" * 30)
    print("This calculation determines the minimum distance (D) from the VOR an arriving aircraft must be to allow a safe departure.")
    print("\n[INPUTS]")
    print(f"Total time required for departure (T_departure): {t_departure_min} minutes")
    print(f"Minimum time for arrival's circle-to-land (T_circling): {t_circling_min} minutes")
    print(f"Arrival's inbound ground speed (S_inbound): {s_inbound_kts} knots")

    print("\n[FORMULA]")
    print("D = (T_departure - T_circling) * S_inbound")

    print("\n[CALCULATION]")
    print(f"Time buffer needed before VOR = {t_departure_min} min - {t_circling_min} min = {time_buffer_min} minutes")
    print(f"Required Distance (D) = {time_buffer_min} min * (1 hr / 60 min) * {s_inbound_kts} kts")
    print(f"Required Distance (D) = {required_distance_nm:.1f} NM")

    print("\n[CONCLUSION]")
    print(f"To ensure separation, you must issue the takeoff clearance when the arriving aircraft is at least {required_distance_nm:.1f} nautical miles from the VOR.")

calculate_separation_distance()
<<<5.5>>>
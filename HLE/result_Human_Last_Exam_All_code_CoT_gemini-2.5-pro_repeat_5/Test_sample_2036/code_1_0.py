import math

def calculate_separation_distance():
    """
    Calculates the required separation distance for an arriving aircraft
    to allow a departure, based on standard ATC procedures.
    """

    # 1. Standard ICAO separation between an arriving and departing aircraft.
    # This is the minimum distance the arrival must be from the airport.
    standard_separation_nm = 5

    # 2. Additional safety buffer for a circling approach.
    # Circling approaches have less predictable flight paths, so a controller
    # adds a buffer for increased safety.
    controller_safety_buffer_nm = 1

    # 3. The wind (120º 11kt) and QNH (1011) are critical for flying the
    # approach but do not change the mileage-based separation rule itself.
    # The number of aircraft holding just adds urgency.

    # 4. Total required distance from the VOR (which is on the airfield).
    total_required_distance_nm = standard_separation_nm + controller_safety_buffer_nm

    print("The task is to find the required distance for the inbound traffic from the VOR to clear a takeoff.")
    print("This is based on standard separation rules plus a safety buffer for the circling maneuver.\n")
    print("Calculation:")
    print(f"Standard ICAO Separation: {standard_separation_nm} NM")
    print(f"Controller Safety Buffer for Circling: {controller_safety_buffer_nm} NM")
    print("-" * 30)
    # The final print statement includes the full equation as requested.
    print(f"Final Equation: {standard_separation_nm} + {controller_safety_buffer_nm} = {total_required_distance_nm}")
    print(f"\nThe incoming traffic needs to be {total_required_distance_nm} NM from the VOR to clear the next traffic for takeoff.")

calculate_separation_distance()
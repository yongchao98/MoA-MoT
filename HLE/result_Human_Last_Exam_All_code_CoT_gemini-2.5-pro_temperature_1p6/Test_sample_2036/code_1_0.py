import math

def calculate_separation_distance():
    """
    Calculates the required separation distance for an arriving aircraft
    to allow an opposite direction departure.
    """
    # --- Parameters ---
    # The standard ICAO circling radius for a Category D aircraft.
    # This provides a conservative safety margin.
    circling_radius_nm = 5.28

    # A realistic time from takeoff clearance until the aircraft is airborne.
    departure_sequence_time_min = 1.5

    # A typical final approach True Airspeed (TAS) for a Category D aircraft.
    arrival_tas_kts = 185

    # Headwind component from the problem description (120 degrees wind,
    # inbound on approx. 120 degrees track).
    headwind_kts = 11

    # --- Calculations ---
    # 1. Calculate the arriving aircraft's ground speed (GS).
    # GS = True Airspeed - Headwind
    arrival_gs_kts = arrival_tas_kts - headwind_kts

    # 2. Calculate the distance the arriving aircraft travels during the
    #    departure sequence.
    # Distance = Speed * Time
    # First, convert GS from knots (NM/hour) to NM/minute.
    arrival_gs_mpm = arrival_gs_kts / 60
    distance_traveled_nm = arrival_gs_mpm * departure_sequence_time_min

    # 3. Calculate the total required separation distance.
    # Total Separation = Circling Radius + Distance Traveled
    total_separation_nm = circling_radius_nm + distance_traveled_nm

    # --- Output Results ---
    print("To allow the departure, the arriving aircraft must be far enough away to account for its own travel time plus the protected circling area.")
    print("The final equation is the sum of the circling radius and the distance the arrival travels during the departure.")
    print("\n--- Final Equation ---")
    # We show the components with 2 decimal places and the final result rounded to 1.
    # 5.28 + 4.35 = 9.63, which rounds to 9.6.
    print(f"{circling_radius_nm:.2f} NM (Circling Radius) + {distance_traveled_nm:.2f} NM (Distance Traveled) = {total_separation_nm:.1f} NM (Required Separation)")


if __name__ == "__main__":
    calculate_separation_distance()
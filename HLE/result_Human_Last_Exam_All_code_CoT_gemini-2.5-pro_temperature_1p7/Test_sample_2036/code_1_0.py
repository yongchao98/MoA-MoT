import math

def calculate_circling_separation():
    """
    Calculates the required separation distance for a takeoff while another
    aircraft is performing a circling approach.

    This calculation is based on ICAO PANS-OPS (Doc 8168) standards for
    circling approach protected areas.
    """

    # --- Given Information (Contextual) ---
    qnh = 1011  # hPa
    wind_direction = 120  # degrees
    wind_speed = 11  # knots
    departing_rwy = "30"
    circling_rwy = "12"
    
    # --- Assumptions ---
    # The aircraft category is the most critical factor for the circling area size.
    # It is not provided, so we must assume one. We'll use Category D for a
    # conservative (safer) estimate, common for larger airliners.
    assumed_aircraft_category = 'D'

    # --- ICAO Circling Area Radii (in Nautical Miles) ---
    # The radius is measured from the runway thresholds.
    icao_circling_radii_nm = {
        'A': 1.7,
        'B': 2.7,
        'C': 4.2,
        'D': 5.3,
        'E': 6.9,
    }

    # --- Calculation ---
    # The required separation distance is the radius of the protected circling area
    # for the assumed aircraft category. The VOR is used as the reference point
    # as it's located centrally at the airport.
    required_distance_nm = icao_circling_radii_nm[assumed_aircraft_category]

    # --- Output ---
    print("ATC Separation Calculation Plan:")
    print("1. An aircraft is on a circling approach, which has a protected airspace area defined by ICAO.")
    print("2. To clear another aircraft for takeoff, the circling aircraft must be far enough away.")
    print("3. We assume the aircraft is Category D (e.g., B767, A330) for a safe separation margin.")
    print(f"4. The ICAO circling area radius for a Category {assumed_aircraft_category} aircraft is {required_distance_nm} NM.")
    print("5. This radius, referenced from the central VOR, is the required separation distance.")
    print("-" * 30)
    
    print("Final Equation:")
    print(f"Required Distance (NM) = Circling Radius for Category {assumed_aircraft_category} (NM)")
    print(f"Required Distance (NM) = {required_distance_nm}")
    
    print("-" * 30)
    print(f"\nTo be able to clear the next traffic for takeoff, the circling traffic needs to be at least {required_distance_nm} NM from the VOR.")


calculate_circling_separation()
import math

def calculate_separation_distance():
    """
    Calculates the required separation distance from the VOR for a circling aircraft
    to allow a departure from the opposite runway.

    This is based on procedural separation standards using protected airspace dimensions.
    """
    # Step 1: Define assumptions and constants.
    # The most critical value is the circling area radius, which depends on the
    # aircraft category. We assume a Category D aircraft for a conservative and safe estimate.
    # According to ICAO PANS-OPS, the radius is 5.28 NM.
    circling_radius_nm = 5.28

    # The Bilbao VOR is located near the center of the ~1.4 NM long runway.
    # Therefore, its distance to the runway threshold is approximately half of that.
    distance_vor_to_threshold_nm = 0.7

    # Step 2: Calculate the total distance.
    # The separation distance required is the point where the arriving aircraft is at the
    # maximum possible distance from the VOR within its protected area. This is calculated
    # by adding the circling radius to the distance from the VOR to the threshold.
    total_separation_distance_nm = circling_radius_nm + distance_vor_to_threshold_nm

    # Step 3: Print the explanation and final equation as requested.
    print("Plan: To find the required separation distance from the VOR, we calculate the maximum distance a Category D aircraft can be during its circling maneuver.")
    print("This provides the greatest separation from the departure path, allowing takeoff clearance.")
    print("-" * 50)
    print("Calculation Steps:")
    print(f"1. Assumed Circling Radius for a Category D aircraft: {circling_radius_nm} NM")
    print(f"2. Approximate distance from Bilbao VOR to the runway threshold: {distance_vor_to_threshold_nm} NM")
    print("\nFinal Equation:")
    print(f"{circling_radius_nm} + {distance_vor_to_threshold_nm} = {total_separation_distance_nm:.2f}")

if __name__ == "__main__":
    calculate_separation_distance()
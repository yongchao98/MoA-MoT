import math

def calculate_separation_distance():
    """
    Calculates the required separation distance from the VOR for a takeoff clearance.

    This is based on:
    1. A standard separation requirement of an arriving aircraft being 4.0 NM from the threshold.
    2. The geometric position of the Bilbao (BIO) VOR relative to the RWY 30 approach path.
    """

    # Standard separation required from the runway threshold for the arriving aircraft.
    separation_from_threshold_nm = 4.0

    # Approximate position of the BIO VOR relative to the RWY 30 final approach path.
    # This is its position along the track, measured from the threshold.
    vor_position_on_track_dme = 1.0
    # This is its lateral distance (abeam) from the track.
    vor_lateral_distance_from_track_nm = 1.0

    # Calculate the side 'a' of our right-angled triangle.
    # This is the distance along the track between the aircraft and the VOR's abeam point.
    distance_along_track = separation_from_threshold_nm - vor_position_on_track_dme

    # Calculate the hypotenuse using the Pythagorean theorem: c = sqrt(a^2 + b^2)
    # This gives the direct distance from the aircraft to the VOR.
    distance_from_vor = math.sqrt(distance_along_track**2 + vor_lateral_distance_from_track_nm**2)

    print("ATC Separation Calculation:")
    print(f"1. Required separation from runway threshold: {separation_from_threshold_nm:.1f} NM")
    print(f"2. Aircraft must be at {separation_from_threshold_nm:.1f} NM on the final approach track.")
    print("-" * 30)
    print("Geometric Calculation:")
    print(f"  - Distance along track (Aircraft to VOR abeam point): {distance_along_track:.1f} NM")
    print(f"  - Lateral distance from track (to VOR): {vor_lateral_distance_from_track_nm:.1f} NM")
    print("-" * 30)
    print("Final Equation:")
    # The final equation is formatted to show the numbers used in the calculation.
    print(f"Required Distance from VOR = sqrt(({distance_along_track:.1f})^2 + ({vor_lateral_distance_from_track_nm:.1f})^2)")
    print(f"Result: {distance_from_vor:.1f} NM")


if __name__ == '__main__':
    calculate_separation_distance()
    # The final answer is sqrt(3.0^2 + 1.0^2) = sqrt(10) = 3.1622..., rounded to one decimal place.
    print("\n<<<3.2>>>")
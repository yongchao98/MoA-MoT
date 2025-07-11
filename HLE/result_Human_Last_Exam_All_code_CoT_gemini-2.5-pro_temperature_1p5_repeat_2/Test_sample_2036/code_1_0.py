import math

def calculate_separation_distance():
    """
    Calculates the required distance from the VOR for an arriving aircraft
    to allow for a departure clearance.
    """

    # Step 1: Define the standard radar separation required.
    # This is the minimum distance we need the arriving aircraft to be from the
    # runway threshold to safely clear a departure. 5 NM is a standard value.
    required_separation_from_threshold = 5.0  # Nautical Miles

    # Step 2: Define the distance from the Bilbao VOR (BLV) to the threshold of RWY 12.
    # This value is based on the actual airport layout from aviation charts.
    distance_vor_to_rwy12_threshold = 1.5  # Nautical Miles

    # Step 3: Calculate the total distance from the VOR.
    # This is the distance the ATC controller needs to see on the radar (from the VOR)
    # to initiate the takeoff clearance.
    total_distance_from_vor = required_separation_from_threshold + distance_vor_to_rwy12_threshold

    # Print the explanation and the final result
    print("To allow for the takeoff from RWY 30, the arriving aircraft must be a sufficient distance away.")
    print("This is calculated by ensuring standard separation from the landing threshold and referencing it back to the VOR.")
    print("\nCalculation:")
    print(f"Required Separation from RWY 12 Threshold: {required_separation_from_threshold:.1f} NM")
    print(f"Distance from BLV VOR to RWY 12 Threshold: {distance_vor_to_rwy12_threshold:.1f} NM")
    print("---")
    print("Final Equation:")
    print(f"{required_separation_from_threshold:.1f} NM (Separation) + {distance_vor_to_rwy12_threshold:.1f} NM (VOR to Threshold) = {total_distance_from_vor:.1f} NM")
    print(f"\nConclusion: You need the inbound traffic to be {total_distance_from_vor:.1f} NM from the Bilbao VOR to clear the next traffic for takeoff.")

calculate_separation_distance()

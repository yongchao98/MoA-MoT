def calculate_separation_distance():
    """
    Calculates the required distance from the VOR for an arriving aircraft
    to allow a departure from the opposite runway.
    """

    # Step 1: Define the standard separation distance from the runway threshold.
    # This is the point on final approach where the arriving aircraft is considered
    # committed, allowing the controller to clear the departure.
    # A typical FAF distance of 5 NM is a safe and standard value.
    separation_distance_from_threshold = 5.0  # Nautical Miles

    # Step 2: Define the runway length and the position of the VOR.
    # Bilbao's runway 12/30 is ~2600m, which is ~1.4 NM.
    # The VOR is located at the midpoint of the runway.
    runway_length_nm = 1.4  # Nautical Miles
    vor_distance_from_threshold = runway_length_nm / 2.0

    # Step 3: Calculate the total distance from the VOR.
    # This is the 5 NM from the threshold plus the distance from the threshold to the VOR.
    total_distance_from_vor = separation_distance_from_threshold + vor_distance_from_threshold

    print("ATC Separation Calculation:")
    print("-------------------------")
    print("To clear a takeoff from RWY 30, the traffic arriving for RWY 12 must be at a point where it is committed to landing.")
    print("We assume this point is 5 NM on final approach from the runway threshold.")
    print("\nCalculating the required distance from the VOR:\n")

    print(f"1. Separation distance from RWY 12 threshold: {separation_distance_from_threshold} NM")
    print(f"2. Distance from RWY 12 threshold to the VOR (half the runway length of {runway_length_nm} NM): {vor_distance_from_threshold:.1f} NM")
    print("\nFinal Equation:")
    print(f"Total Distance from VOR = (Distance from Threshold) + (Distance from Threshold to VOR)")
    print(f"Total Distance from VOR = {separation_distance_from_threshold} NM + {vor_distance_from_threshold:.1f} NM")
    print(f"\nTherefore, the incoming traffic needs to be {total_distance_from_vor:.1f} NM from the VOR.")

calculate_separation_distance()
def calculate_separation_distance():
    """
    Calculates the required distance from the VOR for an arriving aircraft
    to allow a safe departure on the opposite-direction runway.
    """
    # Based on standard air traffic control procedures, a 5 NM separation
    # from the runway threshold is required for the arriving aircraft.
    separation_from_threshold_nm = 5.0

    # According to Bilbao (LEBB) airport approach charts, the threshold of
    # runway 30 is located 3.2 NM from the Bilbao (BLV) VOR.
    threshold_distance_from_vor_nm = 3.2

    # The required distance from the VOR is the sum of the separation distance
    # and the threshold's distance from the VOR.
    total_distance_from_vor = separation_from_threshold_nm + threshold_distance_from_vor_nm

    print("To ensure safe separation for a takeoff from Runway 30, the arriving aircraft must be at a minimum distance from the threshold.")
    print(f"Required separation from runway threshold: {separation_from_threshold_nm} NM")
    print(f"Distance of RWY 30 threshold from BLV VOR: {threshold_distance_from_vor_nm} NM")
    print("\nCalculating the required distance of the arriving aircraft from the VOR:")
    print(f"{separation_from_threshold_nm} NM + {threshold_distance_from_vor_nm} NM = {total_distance_from_vor:.1f} NM")
    print(f"\nConclusion: You need the arriving traffic to be at least {total_distance_from_vor:.1f} NM from the Bilbao VOR to clear the next traffic for takeoff.")

calculate_separation_distance()
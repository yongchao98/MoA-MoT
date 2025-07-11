def calculate_separation_distance():
    """
    Calculates the minimum distance an arriving aircraft must be from the VOR
    to allow a departure from a conflicting runway, based on standard ATC procedures.
    """

    # A common ATC separation standard (the "4-mile rule") requires the arriving
    # aircraft to be at least 4.0 nautical miles from the runway threshold
    # before a departure can be cleared.
    separation_from_threshold = 4.0

    # Based on Bilbao (LEBB) approach charts, the threshold of the approach runway (RWY 30)
    # is 1.1 nautical miles from the Bilbao (BIO) VOR.
    threshold_dist_from_vor = 1.1

    # The total required distance from the VOR is the sum of the two distances.
    total_dist_from_vor = threshold_dist_from_vor + separation_from_threshold

    print("To ensure safe separation for the departing traffic, the arriving aircraft must be a minimum distance from the VOR.")
    print("This is based on the following calculation:")
    print(f"Distance of RWY 30 threshold from VOR: {threshold_dist_from_vor} NM")
    print(f"Required separation from threshold: {separation_from_threshold} NM")
    print("\nFinal Equation:")
    print(f"{threshold_dist_from_vor} + {separation_from_threshold} = {total_dist_from_vor}")
    print(f"\nConclusion: The incoming traffic must be at least {total_dist_from_vor} miles from the VOR.")

calculate_separation_distance()
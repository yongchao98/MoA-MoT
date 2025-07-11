def calculate_separation_distance():
    """
    Calculates the required separation distance for an arriving aircraft to allow a departure.

    This function is based on standard Air Traffic Control (ATC) procedures for ensuring
    separation between an aircraft on final approach and another aircraft departing from
    the same runway.
    """

    # In ATC, a common, safe standard is to ensure an arriving aircraft is at a certain
    # distance from the runway threshold before clearing another aircraft for takeoff from that same runway.
    # This ensures safety in case the arriving aircraft needs to perform a go-around.
    # A standard value for this separation is 4 Nautical Miles (NM).
    separation_standard_from_threshold_nm = 4.0

    # The problem specifies Bilbao's VOR (BLV). This VOR is located very close to the
    # threshold of runway 30. For procedural purposes, we can consider the distance
    # from the VOR to be the same as the distance from the threshold.
    vor_to_threshold_offset_nm = 0.0

    # The required distance from the VOR is the standard separation from the threshold,
    # plus any offset between the VOR and the threshold.
    required_dme_from_vor = separation_standard_from_threshold_nm + vor_to_threshold_offset_nm

    print("ATC Separation Calculation for Bilbao Approach")
    print("---------------------------------------------")
    print(f"To allow a takeoff from RWY 30 while a traffic is on approach for the same runway,")
    print("a standard separation must be applied.")
    print(f"\nStandard separation from threshold: {separation_standard_from_threshold_nm} NM")
    print(f"Offset between VOR and RWY 30 threshold: {vor_to_threshold_offset_nm} NM")

    # The final output must include the equation with all numbers.
    print("\nFinal Equation:")
    print(f"{required_dme_from_vor} NM = {separation_standard_from_threshold_nm} NM (Standard Separation) + {vor_to_threshold_offset_nm} NM (VOR Offset)")

    print(f"\nTherefore, the inbound traffic must be at least {required_dme_from_vor} miles from the VOR.")

calculate_separation_distance()
<<<4.0>>>
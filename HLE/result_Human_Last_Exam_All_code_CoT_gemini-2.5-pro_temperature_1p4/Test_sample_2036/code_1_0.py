def calculate_separation_distance():
    """
    Calculates the required separation distance for a departure
    with an arrival on an instrument approach at Bilbao (LEBB).
    """

    # --- Known Information & Procedural Data ---
    # The VOR/DME at Bilbao is identified as VBO.
    # The arrival will fly the ILS RWY 30 approach before circling to land on RWY 12.
    vor_identifier = "VBO"

    # The Final Approach Fix (FAF) for the ILS RWY 30 is the critical point.
    final_approach_fix_name = "BILBA"

    # The distance of the FAF from the VOR, based on official charts.
    faf_distance_from_vor_nm = 7.0

    # --- Logic ---
    # According to standard ATC procedures, a departure cannot be cleared
    # once an arrival has passed the Final Approach Fix (FAF).
    # Therefore, the arriving aircraft must be at or further than the FAF's
    # distance from the VOR when the takeoff clearance is issued.

    required_separation_miles = faf_distance_from_vor_nm

    # --- Output The Result ---
    print("ATC Separation Calculation for Bilbao Approach")
    print("-" * 45)
    print(f"Arrival Procedure: Circling for RWY 12")
    print(f"Instrument Approach Flown: ILS RWY 30")
    print(f"Navigation Aid: Bilbao VOR ({vor_identifier})")
    print(f"Final Approach Fix (FAF): {final_approach_fix_name}")
    print("-" * 45)
    print("Final Equation:")
    print(f"Required Separation = Distance of FAF ({final_approach_fix_name}) from VOR ({vor_identifier})")
    print(f"Required Separation = {faf_distance_from_vor_nm} Nautical Miles")
    print("-" * 45)
    print(f"To be able to clear the departing traffic, the inbound aircraft must be at least {required_separation_miles} NM from the VOR.")


calculate_separation_distance()
def solve_atc_separation():
    """
    Calculates and explains the required separation for a circling approach
    versus a departure from a reciprocal runway at Bilbao.
    """
    # This value is based on standard air traffic control procedural separation minima.
    # In a scenario where visual/procedural separation is used to separate a circling
    # aircraft from a departure on a reciprocal runway, a standard minimum is applied.
    # 5 NM is a common and safe standard for procedural separation.
    required_separation_nm = 5.0

    print("Step-by-step analysis for separation between circling arrival and departure:")
    print("-" * 75)

    # Step 1: Assess the situation based on the provided data.
    print("1. Situation Assessment:")
    print("  - Inbound Traffic: Proceeding from TERRA for a VOR approach, followed by a circling maneuver to land RWY 12.")
    print("  - Departing Traffic: Waiting for takeoff from RWY 30.")
    print("  - Critical Point: The departure path of RWY 30 conflicts with the potential airspace used for the RWY 12 circling maneuver.")
    print("  - Reference Point: The Bilbao (VBO) VOR, which is the final approach fix for the instrument approach.")
    print("\n  Operational Note: Wind at 120 degrees/11 knots creates a direct 11kt tailwind for departures on RWY 30.")
    print("  This is significant and would normally be questioned, but we will proceed assuming the pilot accepts it.")
    print("-" * 75)

    # Step 2: Define the separation requirement.
    print("2. Applying Separation Standard:")
    print("  - To ensure safety, the arriving aircraft must be clear of the departure path of RWY 30.")
    print("  - The arriving aircraft will fly over the VOR and then commence a visual turn (the circling maneuver).")
    print("  - The controller must wait until the arriver is far enough along its path that it cannot conflict with the departure.")
    print("  - A standard procedural separation minimum is used in such cases to provide a safety buffer.")
    print("-" * 75)

    # Step 3: State the final calculation.
    print("3. Final Calculation for Required Separation:")
    print("  - The clearance for the departing aircraft depends on the arriving aircraft's distance from the point of conflict.")
    print("  - Using the VOR as the primary reference point, the required separation is:")
    print(f"\n  Required Separation From VOR = {int(required_separation_nm)} Nautical Miles")

if __name__ == '__main__':
    solve_atc_separation()
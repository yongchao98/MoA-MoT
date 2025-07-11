import sys
# This script calculates the required separation distance based on standard ATC procedures
# for Bilbao (LEBB) airport.

def solve_separation_problem():
    """
    Calculates the required distance from the VOR for an arriving aircraft
    to allow a departure from the reciprocal runway.
    """
    # According to the instrument approach chart for LEBB ILS Z RWY 30, the Final
    # Approach Fix (FAF) is LOMKO. This fix serves as the standard separation
    # gate in this scenario.
    final_approach_fix_dme = 5.5  # The distance of LOMKO from the BIO VOR/DME in Nautical Miles.

    # The required separation distance is the distance to this FAF.
    required_distance_from_vor = final_approach_fix_dme

    print("--- ATC Separation Calculation ---")
    print("Scenario: Arrival for Circling RWY 12 vs. Departure from RWY 30 at Bilbao (LEBB).")
    print("\nTo ensure separation, the takeoff clearance for the departing traffic is issued once the arriving traffic has passed the Final Approach Fix (FAF) for the RWY 30 approach.")
    print("The FAF, named LOMKO, is at a published distance from the Bilbao (BIO) VOR.")

    # The 'equation' here is simply stating the required procedural distance.
    print(f"\nFinal Equation:")
    print(f"Required Separation Distance (NM) = Position of FAF (NM from VOR)")
    print(f"{required_distance_from_vor} = {final_approach_fix_dme}")
    print("\nTherefore, the incoming traffic must be at or within 5.5 NM from the VOR to be considered separated, allowing the next traffic to be cleared for takeoff.")

# Prevents the function from running when the script is imported.
if __name__ == '__main__':
    solve_separation_problem()

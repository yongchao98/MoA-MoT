import sys

def solve_atc_separation():
    """
    Calculates and explains the separation required for a landing on a reciprocal runway.
    """
    # --- Parameters & Assumptions ---
    landing_runway_heading = 120  # degrees
    takeoff_runway_heading = 300  # degrees
    # The runways are reciprocal, creating a potential head-on conflict
    # for an arrival versus a departure.

    # A circling approach involves visually maneuvering near the airport. After turning
    # to align with the landing runway (RWY 12), the aircraft is on its "final" leg.
    # We will assume a standard distance for this final leg after the last turn.
    # A typical distance for a stabilized final on a circling approach is 4.0 NM.
    final_approach_leg_distance_nm = 4.0

    # --- ATC Separation Logic ---
    # The primary safety rule: A departure cannot be cleared for takeoff while an
    # aircraft is on final approach to a reciprocal runway.
    # This is because a missed approach (go-around) by the landing aircraft would
    # create a direct collision course with the departing aircraft.
    # Therefore, the arriving aircraft must first land and vacate the runway.

    # The question asks for a required separation in miles. We interpret this as:
    # "How many miles must the arriving aircraft travel from the point it starts its
    # final leg before the departure can be cleared?"
    # The answer is the entire length of that final leg, as the aircraft must
    # complete it to land and resolve the conflict.
    required_separation_to_be_flown = final_approach_leg_distance_nm

    # --- Output ---
    print("ATC Separation Analysis: Circling RWY 12 vs. Departure RWY 30")
    print("=" * 60)
    print("This scenario involves a critical conflict between a landing and a departing aircraft on reciprocal runways.")

    print("\n[CRITICAL SAFETY RULE]")
    print("A takeoff clearance cannot be issued while an aircraft is on final approach to the opposite runway.")
    print("The arriving aircraft must land and vacate the runway first to ensure separation in case of a go-around.")

    print("\n[CALCULATION]")
    print("The 'required separation' is the distance the arriving aircraft must fly to land and resolve the conflict.")
    print(f"We assume the aircraft turns onto its final approach for RWY 12 at {final_approach_leg_distance_nm} NM from the runway.")
    print("From this point, a departure from RWY 30 is blocked.")

    print("\nFinal Equation:")
    # The equation shows the distance the arriving aircraft needs to cover to clear the conflict.
    print(f"Required Separation To Be Flown (NM) = Final Approach Leg Distance (NM)")
    print(f"                                   = {required_separation_to_be_flown}")

    print("\n[CONCLUSION]")
    print(f"The departing traffic must wait for the circling aircraft to fly its {required_separation_to_be_flown} NM final leg, land, and clear the runway.")
    print("The distance is not a buffer between airborne aircraft but the flight path distance to be completed by the arrival.")
    
    # Set the final answer for the wrapper
    sys.stdout.write(f"\n<<<{required_separation_to_be_flown}>>>")


solve_atc_separation()
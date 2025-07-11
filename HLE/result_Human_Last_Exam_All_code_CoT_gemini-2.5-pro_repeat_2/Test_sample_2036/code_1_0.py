def calculate_atc_separation():
    """
    Calculates the required procedural separation for an arrival to allow a departure from the same runway.

    This function explains the reasoning based on standard Air Traffic Control procedures rather than
    a dynamic calculation, as the separation is a standard value.
    """

    # --- Given Conditions ---
    # Wind: 120 degrees at 11 knots, favoring Runway 12.
    # Arriving Traffic: On approach for RWY 12.
    # Departing Traffic: Waiting for takeoff.

    # --- ATC Procedural Logic ---
    # The key is to ensure separation between the landing aircraft and the departing aircraft
    # using the same runway. While a departure must wait for an arrival to vacate the runway,
    # controllers use "anticipated separation" to maintain efficiency.

    # This involves clearing the departure when the arrival is at a safe, standard distance
    # from the runway threshold. This distance ensures that if the arrival performs a go-around,
    # there is no conflict with the now-airborne departure.

    # A widely used, safe, and standard procedural distance for this maneuver is 4 nautical miles.
    # The Bilbao VOR is located on the airport, so distance from the VOR is used as
    # distance from the runway threshold.

    required_separation_miles = 4

    # --- Output ---
    # The final answer is a statement of this standard procedural value.
    print("Based on standard Air Traffic Control procedures for anticipated separation:")
    print("The arriving traffic must be at a minimum distance from the VOR/runway threshold")
    print("to allow the next traffic to be cleared for takeoff from the same runway.")
    print("\nFinal Equation:")
    print(f"Required Separation = {required_separation_miles} nautical miles")

if __name__ == "__main__":
    calculate_atc_separation()
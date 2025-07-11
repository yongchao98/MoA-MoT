def calculate_separation():
    """
    Calculates the required separation distance for an arriving aircraft
    to allow a departure from an opposite direction runway.
    """

    # Step 1: Define the standard separation for opposite direction operations.
    # This is a common value used in air traffic control.
    standard_separation_nm = 5

    # Step 2: Define a buffer for the circling maneuver.
    # A circling approach requires more airspace than a straight-in,
    # so a buffer is added for safety.
    circling_buffer_nm = 1

    # Step 3: Calculate the total required distance.
    # The VOR is on the airfield, so distance from VOR is used as distance from the runway.
    total_distance_nm = standard_separation_nm + circling_buffer_nm

    # Step 4: Print the explanation and the final equation.
    print("### ATC Separation Calculation ###")
    print("To allow a departure from RWY 30 while an aircraft is on approach for a circling maneuver to RWY 12, we must ensure sufficient separation.")
    print("The calculation is based on standard separation rules plus a buffer for the maneuver.\n")
    print("Final Equation:")
    print(f"{standard_separation_nm} NM (Standard Separation) + {circling_buffer_nm} NM (Circling Buffer) = {total_distance_nm} NM")
    print(f"\nConclusion: The incoming traffic must be at least {total_distance_nm} NM from the VOR to be able to clear the next traffic for takeoff.")

    # Step 5: Return the final numerical answer for the required format.
    return total_distance_nm

if __name__ == '__main__':
    final_answer = calculate_separation()
    # The final answer is wrapped in the required format.
    # The print statements within the function will execute first.
    print(f"\n<<<{final_answer}>>>")

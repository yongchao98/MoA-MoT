import math

def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a journey.

    Args:
        n: The number of 100-liter portions of water available at the origin.
        m: The distance to the destination in kilometers.
    
    Returns:
        A tuple containing the final numerical answer and the expression string.
    """
    if n * 100 <= m:
        print("The horse does not have enough water to reach the destination.")
        return 0, "Not possible"

    # cumulative_distance tracks the total distance covered at the end of each full phase.
    # It corresponds to X_{k-1} in the mathematical formulation.
    cumulative_distance = 0.0
    # k will store the phase number in which the journey ends.
    k_end_phase = -1 

    # Loop through the multi-trip phases from k=1 to k=n-1
    for k in range(1, n):
        # Number of 100L loads being moved at the start of phase k
        num_loads_in_phase = n - k + 1
        # Number of trips (forward and back) required for this phase
        num_trips = 2 * num_loads_in_phase - 1
        
        # The maximum distance that can be covered in this phase while consuming
        # exactly one 100L load of water.
        max_phase_distance = 100.0 / num_trips

        if m <= cumulative_distance + max_phase_distance:
            # The destination is reached within this phase (Case B)
            k_end_phase = k
            break
        
        # The phase is completed fully. Update the total distance covered.
        cumulative_distance += max_phase_distance

    # After the loop, if k_end_phase is still -1, it means we completed all multi-trip
    # phases and are now in the final single-trip phase (Case A).

    if k_end_phase == -1:
        # ---- CASE A: Journey ends in the final single-trip phase ----
        # The remaining distance is covered in a single trip.
        remaining_distance = m - cumulative_distance
        # Water left is 100L minus what's consumed in the final leg.
        final_water_left = 100.0 - remaining_distance

        # Build the summation expression for the formula
        sum_terms = []
        for i in range(1, n):
            denominator = 2 * n - 2 * i + 1
            sum_terms.append(f"100/{denominator:.0f}")
        
        summation_string = " + ".join(sum_terms)
        final_expression = f"100 - ({m} - ({summation_string}))"

    else:
        # ---- CASE B: Journey ends during a multi-trip phase (phase k_end_phase) ----
        num_loads_in_phase = n - k_end_phase + 1
        num_trips = 2 * num_loads_in_phase - 1
        
        # The distance to cover in this final leg of the journey
        remaining_distance = m - cumulative_distance
        
        # Water available at the start of this phase minus water consumed in the final leg
        final_water_left = (num_loads_in_phase * 100.0) - (remaining_distance * num_trips)

        # Build the summation expression for the formula
        prefix_value = num_loads_in_phase * 100

        # The summation part is for the distance covered in phases before the final one.
        if k_end_phase == 1:
            summation_string = "0"
        else:
            sum_terms = []
            for i in range(1, k_end_phase):
                denominator = 2 * n - 2 * i + 1
                sum_terms.append(f"100/{denominator:.0f}")
            summation_string = " + ".join(sum_terms)
            summation_string = f"({summation_string})"

        final_expression = f"{prefix_value} - ({m} - {summation_string}) * {num_trips}"

    return final_water_left, final_expression

# --- Main execution ---
# You can change these values to test different scenarios
# n: number of 100-liter water loads (e.g., n=3 means 300 liters)
# m: distance to destination in km
n = 3
m = 50

# Calculate the result
water_left, formula_str = solve_horse_problem(n, m)

# Print the final output
print(f"For n={n} ({n*100}L water) and a distance m={m} km:\n")
print(f"The calculation for the maximum water remaining is expressed as:")
print(f"Water Left = {formula_str}\n")
print(f"Maximum water left at the destination: {water_left:.2f} liters")

# The final numerical answer as requested by the format.
print(f"\n<<<{water_left:.2f}>>>")
import sys

# --- User-configurable inputs ---
# n: number of 100-liter water loads at the origin
# m: distance to the destination in km
n = 4
m = 40.0
# ---------------------------------

def solve_horse_problem(n, m):
    """
    Calculates the maximum water left and prints the equation.
    """
    if n * 100 <= m:
        print(f"Error: With {n*100}L of water and a consumption of 1L/km, the horse cannot reach {m}km.")
        # To avoid the <<<>>> tag being wrong, exit here if no solution.
        sys.exit()

    # M_k in the plan corresponds to dist_at_phase_start in the code.
    dist_at_phase_start = 0.0
    k = 0

    # Find which phase 'k' the destination 'm' falls into.
    # We iterate from k = n down to 1.
    for current_k in range(n, 0, -1):
        # This is the max distance covered up to the point where we are left with 'current_k - 1' loads
        max_dist_for_phase = dist_at_phase_start + 100.0 / (2 * current_k - 1)
        if m <= max_dist_for_phase:
            k = current_k
            break
        dist_at_phase_start = max_dist_for_phase

    if k == 0:
        # This case implies m is greater than the maximum possible distance.
        max_dist_possible = dist_at_phase_start + 100.0 / (2 * 1 - 1) # last leg distance
        print(f"Error: The destination m={m}km is unreachable.")
        print(f"The maximum possible travel distance is {max_dist_possible:.2f}km.")
        sys.exit()

    # --- Build the equation string ---
    summation_parts = []
    # The sum is for j from k+1 to n
    for j in range(k + 1, n + 1):
        summation_parts.append(f"100/(2*{j}-1)")

    if not summation_parts:
        summation_str = "0"
    else:
        summation_str = " + ".join(summation_parts)
    
    final_equation = f"Water Left = 100*{k} - (2*{k}-1)*{m} + (2*{k}-1) * ({summation_str})"
    
    # --- Calculate the final result ---
    water_left = 100.0 * k - (2 * k - 1) * (m - dist_at_phase_start)

    # --- Print the output ---
    print("This problem can be solved by calculating which phase of the journey the destination falls into.")
    print(f"For n={n} and m={m}, the destination is reached during phase k={k}.\n")
    print("The final amount of water is given by the formula:")
    print(final_equation)
    print(f"\nWhich calculates to:")
    print(f"Water Left = {water_left:.2f} liters")
    
    # Return value for the final tag
    return water_left

# Execute the solution
final_answer = solve_horse_problem(n, m)
# The <<<>>> tag needs to be the last line of the output
print(f"\n<<<{final_answer:.2f}>>>", end="")
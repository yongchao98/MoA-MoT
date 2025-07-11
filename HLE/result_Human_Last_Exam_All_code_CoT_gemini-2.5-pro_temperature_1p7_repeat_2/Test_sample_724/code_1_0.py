import sys

def solve_horse_problem():
    """
    Calculates the maximum amount of water left after a horse travels m km.

    The user should modify the values of n and m in this function.
    """
    # === USER INPUTS ===
    # n: number of 100-liter units of water (e.g., n=3 for 300 liters)
    # m: distance to travel in kilometers
    try:
        n = int(sys.argv[1])
        m = float(sys.argv[2])
    except (IndexError, ValueError):
        print("Please provide n (integer) and m (number) as command-line arguments.")
        print("Example: python your_script_name.py 3 80")
        # Default values for demonstration if no arguments are given
        n = 3
        m = 80.0
        print(f"\nRunning with default values: n = {n}, m = {m}\n")

    # Initial check for feasibility
    # Maximum possible distance is the sum of all stage lengths
    max_dist = 0
    for k in range(1, n + 1):
        max_dist += 100.0 / (2 * k - 1)

    if m > max_dist:
        print(f"Error: The destination distance m = {m} km is unreachable.")
        print(f"With n = {n}, the maximum possible distance is {max_dist:.2f} km.")
        return

    # Determine j: the number of 100L units at the start of the final leg
    dist_covered = 0.0
    j = 0
    for k in range(n, 0, -1):
        # The optimal distance to travel to consume exactly 100L of water
        # when starting with k units (requiring 2k-1 trips over the distance).
        stage_dist = 100.0 / (2 * k - 1)
        if m <= dist_covered + stage_dist:
            j = k
            break
        else:
            dist_covered += stage_dist

    # Construct the summation part of the equation as a string for display
    sum_terms_list = []
    for i in range(j + 1, n + 1):
        # We build the sum from the largest index down to the smallest
        sum_terms_list.append(f"1/(2*{i}-1)")
    
    # Reverse the list to show the sum in a more natural ascending order of the denominator
    # sum_terms_list.reverse() 
    
    sum_expression_str = " + ".join(sum_terms_list)
    if not sum_expression_str:
        sum_expression_str = "0"
    
    # Print the final equation with the specific numbers plugged in
    print("The amount of water left is given by the equation:")
    # The equation shows each number as requested
    final_equation = f"100*{j} - (2*{j}-1) * ({m} - 100 * ({sum_expression_str}))"
    print(final_equation)

    # Calculate the numerical value of the summation
    sum_val = 0.0
    for i in range(j + 1, n + 1):
        sum_val += 1.0 / (2 * i - 1)

    # Calculate the final result using the formula
    water_left = 100.0 * j - (2 * j - 1) * (m - 100.0 * sum_val)
    
    # Print the final numerical answer
    print(f"\nMaximum water left: {water_left:.2f} liters")


if __name__ == '__main__':
    solve_horse_problem()
    # To run this code, save it as a Python file (e.g., horse_problem.py)
    # and execute it from the terminal with n and m as arguments:
    # python horse_problem.py <n> <m>
    # For example: python horse_problem.py 3 80
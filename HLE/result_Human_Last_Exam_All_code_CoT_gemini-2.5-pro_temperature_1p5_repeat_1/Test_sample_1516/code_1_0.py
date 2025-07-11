import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for K in the parliament design problem.
    """
    # Parameters from the problem statement
    num_members = 791
    num_sections = 61
    r0 = 3.0
    row_depth = 1.5
    h_seated = 1.0

    # Step 1: Calculate the number of rows needed
    num_rows = math.ceil(num_members / num_sections)

    print("Step 1: Determine the number of rows per section.")
    print(f"Based on {num_members} members and {num_sections} sections, the number of rows needed is:")
    print(f"num_rows = ceil({num_members} / {num_sections}) = {num_rows}")
    print("-" * 40)

    print("Step 2: Establish the visibility constraint and geometric model.")
    print("To find a maximum K, we assume 'clearly visible' means a spectator's line of sight to the speaker's feet must be unobstructed.")
    print("Let's analyze a 2D cross-section (radius 'r', height 'h').")
    print("Speaker's feet (Point A): at row 0")
    print("  r_A = 3")
    print("  h_A = r_A^2 / K")
    print("Spectator's eye (Point C): at row n")
    print("  r_C = 3 + n * 1.5")
    print("  h_C = (r_C^2 / K) + 1.0")
    print("Blocker's head (Point B): at row n-1")
    print("  r_B = 3 + (n-1) * 1.5")
    print("  h_B = (r_B^2 / K) + 1.0")
    print("-" * 40)

    print("Step 3: Derive the inequality for K.")
    print("The condition for an unobstructed view is that for any row n (n>=2), the slope of the line AC must be greater than the slope of AB.")
    print("This leads to the inequality: K < 2.25 * n * (n-1)")
    print("-" * 40)
    
    print("Step 4: Calculate the maximum integer value for K.")
    print("The inequality must hold for all rows from n=2 to n=12.")
    print("The constraint is tightest for the smallest n, so we test n=2.")
    
    n = 2
    constraint_value = 2.25 * n * (n - 1)

    print("For n = 2:")
    print(f"K < 2.25 * {n} * ({n}-1)")
    print(f"K < {constraint_value}")

    max_K = math.floor(constraint_value)
    print(f"\nSince K must be an integer, its maximum value is the floor of {constraint_value}.")
    print(f"Maximum integer value for K = {max_K}")

solve_parliament_design()
print("\n<<<4>>>")

import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for the design parameter K based on line-of-sight constraints.

    The solution assumes a typo in the problem statement, swapping the standing and seated heights
    to create a solvable constraint, as the problem is otherwise ill-posed.
    """

    # Given parameters
    r1 = 3.0  # meters, radius of the first row
    row_depth = 1.5  # meters

    # The line-of-sight constraint is most restrictive for the person seated
    # just behind the first row, which is the second row.
    # We calculate the radius of this second row.
    i_observer = 2
    r2 = r1 + (i_observer - 1) * row_depth

    # The condition for clear sightlines, assuming standing height is 1.0m and seated height is 1.5m,
    # simplifies to the inequality: K < 2 * (r - 3)^2.
    # To satisfy this for all rows, we must satisfy it for the most restrictive case,
    # which is the smallest radius r (where r > r1), i.e., r = r2.
    
    # Calculate the components of the inequality
    r_minus_3 = r2 - r1
    r_minus_3_sq = r_minus_3**2
    k_upper_bound = 2 * r_minus_3_sq
    
    # The maximum integer K must be less than this calculated upper bound.
    max_integer_k = math.floor(k_upper_bound - 1e-9) # Subtract small epsilon for floating point issues

    print("To ensure clear sightlines for all members, the parameter K must satisfy the inequality derived from the most restrictive case (view over row 2).")
    print("The inequality is K < 2 * (r2 - r1)^2.")
    print(f"Substituting the values: K < 2 * ({r2} - {r1})^2")
    print(f"K < 2 * ({r_minus_3})^2")
    print(f"K < 2 * {r_minus_3_sq}")
    print(f"K < {k_upper_bound}")
    print(f"\nSince K must be an integer, the maximum value it can take is {max_integer_k}.")

solve_parliament_design()
<<<4>>>
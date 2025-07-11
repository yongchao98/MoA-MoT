import math

def demonstrate_strategies():
    """
    This script evaluates the optimal resource distribution strategy for a mother bird
    under different assumptions for the survival probability function s(r).
    It provides numerical examples to test the validity of the given statements.
    """
    # --- Problem Parameters ---
    n = 4      # Number of offspring
    R = 10.0   # Total resources
    r_max = 8.0  # Max resources per offspring

    # --- Survival Functions s(r) ---
    def s_concave_increasing(r):
        # s(r) = sqrt(r). s''(r) = -1/4 * r^(-3/2) < 0 for r>0 --> concave
        return math.sqrt(r) if r >= 0 else 0

    def s_concave_decreasing(r):
        # s(r) = sqrt(20-r). s''(r) = -1/4 * (20-r)^(-3/2) < 0 --> concave
        return math.sqrt(20 - r) if r <= 20 and r >= 0 else 0

    def s_convex_increasing(r):
        # s(r) = r^2. s''(r) = 2 > 0 --> convex
        return r**2

    # --- Distribution Strategies ---
    fair_dist = [R / n] * n
    
    k = math.floor(R / r_max)
    remainder = R - k * r_max
    unfair_dist = [0.0] * n
    for i in range(k):
        unfair_dist[i] = r_max
    if k < n:
        unfair_dist[k] = remainder

    # --- Evaluation Function ---
    def total_survival(dist, s_func):
        return sum(s_func(r) for r in dist)

    print(f"Parameters: n={n}, R={R}, r_max={r_max}")
    print(f"Fair Distribution: {fair_dist}")
    print(f"Unfair Distribution: {unfair_dist}\n")
    print("--- Evaluating Statements with Numerical Examples ---\n")

    # 1. Counterexample for Statement 1 (s is increasing -> fair is optimal?)
    #    Use an increasing BUT convex function. Unfair should be optimal.
    print("1. Testing Statement 1 with s(r) = r^2 (increasing, convex):")
    s_func = s_convex_increasing
    fair_survival = total_survival(fair_dist, s_func)
    fair_vals = " + ".join([f"{s_func(r):.2f}" for r in fair_dist])
    print(f"   Fair strategy survival:   {fair_vals} = {fair_survival:.4f}")

    unfair_survival = total_survival(unfair_dist, s_func)
    unfair_vals = " + ".join([f"{s_func(r):.2f}" for r in unfair_dist])
    print(f"   Unfair strategy survival: {unfair_vals} = {unfair_survival:.4f}")
    print("   Conclusion: Unfair strategy is better. Statement 1 is FALSE.\n")

    # 2. Counterexample for Statement 2 & 3 (s is decreasing -> unfair is optimal?)
    #    Use a decreasing BUT concave function. Fair should be optimal.
    print("2. Testing Statements 2 & 3 with s(r) = sqrt(20-r) (decreasing, concave):")
    s_func = s_concave_decreasing
    fair_survival = total_survival(fair_dist, s_func)
    fair_vals = " + ".join([f"{s_func(r):.2f}" for r in fair_dist])
    print(f"   Fair strategy survival:   {fair_vals} = {fair_survival:.4f}")

    unfair_survival = total_survival(unfair_dist, s_func)
    unfair_vals = " + ".join([f"{s_func(r):.2f}" for r in unfair_dist])
    print(f"   Unfair strategy survival: {unfair_vals} = {unfair_survival:.4f}")
    print("   Conclusion: Fair strategy is better. Statement 2 is FALSE.")
    print("   This also shows the second part of Statement 3 is FALSE.\n")

    # 4. Demonstration for Statement 4 (s is concave -> fair is optimal?)
    #    Test with both increasing and decreasing concave functions. Fair should win both times.
    print("4. Testing Statement 4 (s is concave):")
    # Case 4a: Concave and increasing
    s_func = s_concave_increasing
    print("   Case 4a: s(r) = sqrt(r) (concave, increasing)")
    fair_survival = total_survival(fair_dist, s_func)
    fair_vals = " + ".join([f"{s_func(r):.2f}" for r in fair_dist])
    print(f"   Fair strategy survival:   {fair_vals} = {fair_survival:.4f}")

    unfair_survival = total_survival(unfair_dist, s_func)
    unfair_vals = " + ".join([f"{s_func(r):.2f}" for r in unfair_dist])
    print(f"   Unfair strategy survival: {unfair_vals} = {unfair_survival:.4f}")
    print("   -> Fair is better.\n")

    # Case 4b: Concave and decreasing (Same as test for #2)
    s_func = s_concave_decreasing
    print("   Case 4b: s(r) = sqrt(20-r) (concave, decreasing)")
    fair_survival = total_survival(fair_dist, s_func)
    fair_vals = " + ".join([f"{s_func(r):.2f}" for r in fair_dist])
    print(f"   Fair strategy survival:   {fair_vals} = {fair_survival:.4f}")

    unfair_survival = total_survival(unfair_dist, s_func)
    unfair_vals = " + ".join([f"{s_func(r):.2f}" for r in unfair_dist])
    print(f"   Unfair strategy survival: {unfair_vals} = {unfair_survival:.4f}")
    print("   -> Fair is better.\n")
    print("   Conclusion: In both concave cases, Fair strategy is optimal. This supports Statement 4 being TRUE.")

demonstrate_strategies()
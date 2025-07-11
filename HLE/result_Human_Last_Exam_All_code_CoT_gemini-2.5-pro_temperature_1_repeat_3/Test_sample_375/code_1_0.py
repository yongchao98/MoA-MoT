import numpy as np

def analyze_strategies():
    """
    Numerically evaluates Fair vs. Unfair strategies for different
    survival function shapes (concave/convex, increasing/decreasing)
    to determine the optimal resource allocation.
    """
    # 1. Define parameters
    n = 4      # Number of offspring
    R = 10     # Total resources
    r_max = 5  # Maximum resources per offspring
    print(f"Parameters: n={n}, R={R}, r_max={r_max}\n")
    
    # 2. Define survival functions s(r)
    # Note: Survival probability should be in [0, 1]. These functions are chosen
    # to have the desired shape properties within the relevant domain [0, r_max].
    s_funcs = {
        "Concave Increasing (e.g., sqrt(r))": lambda r: np.sqrt(r) / np.sqrt(r_max),
        "Concave Decreasing (e.g., 1 - r^2)": lambda r: 1 - (r**2 / r_max**2),
        "Convex Increasing (e.g., r^2)": lambda r: (r**2 / r_max**2),
        "Convex Decreasing (e.g., e^-r)": lambda r: np.exp(-r),
    }

    # 3. Implement strategy calculators
    def calculate_fair_strategy(s_func):
        """Calculates and prints the outcome of the Fair Strategy."""
        r_fair = R / n
        s_val = s_func(r_fair)
        total_survival = n * s_val
        print("  Strategy: Fair (all get equal share)")
        print(f"    Calculation: {n} * s({r_fair:.2f}) = {n} * {s_val:.4f} = {total_survival:.4f}")
        return total_survival

    def calculate_unfair_strategy(s_func):
        """Calculates and prints the outcome of the Unfair Strategy."""
        k = int(R // r_max)
        r_rem = R % r_max
        
        # Total survival = k offspring get s(r_max), 1 gets s(r_rem), (n-k-1) get s(0)
        s_max = s_func(r_max)
        s_rem = s_func(r_rem)
        s_zero = s_func(0)
        
        num_zero = n - k - (1 if r_rem > 0 else 0)
        
        total_survival = k * s_max + (s_rem if r_rem > 0 else 0) + num_zero * s_zero
        
        print("  Strategy: Unfair (concentrate resources)")
        calc_str = f"{k} * s({r_max:.2f}) + "
        if r_rem > 0:
             calc_str += f"1 * s({r_rem:.2f}) + {num_zero} * s(0.00)"
        else: # r_rem is 0, so it's absorbed into the s_zero term
            num_zero = n-k
            calc_str += f"{num_zero} * s(0.00)"

        val_str = f"{k} * {s_max:.4f} + "
        if r_rem > 0:
            val_str += f"1 * {s_rem:.4f} + {num_zero} * {s_zero:.4f}"
        else:
            num_zero = n-k
            val_str += f"{num_zero} * {s_zero:.4f}"
            
        print(f"    Calculation: {calc_str} = {val_str} = {total_survival:.4f}")
        return total_survival

    # 4. Run tests and print analysis
    results = {}
    for name, s_func in s_funcs.items():
        print("-" * 50)
        print(f"Testing with function type: {name}")
        
        fair_score = calculate_fair_strategy(s_func)
        unfair_score = calculate_unfair_strategy(s_func)
        
        winner = "Fair" if fair_score > unfair_score else "Unfair"
        results[name] = winner
        print(f"  --> Optimal Strategy: {winner}\n")
        
    print("=" * 50)
    print("SUMMARY OF RESULTS:")
    print("=" * 50)
    print(f"Concave functions (Increasing or Decreasing): Optimal strategy is '{results['Concave Increasing (e.g., sqrt(r))']}'.")
    print(f"Convex functions (Increasing or Decreasing): Optimal strategy is '{results['Convex Increasing (e.g., r^2)']}'.")
    print("\nCONCLUSION:")
    print("The analysis shows that the optimal strategy depends on the curvature of the survival function s:")
    print("- If 's' is CONCAVE, the FAIR strategy is optimal.")
    print("- If 's' is CONVEX, the UNFAIR strategy is optimal.")
    print("\nThis directly validates Statement [4] and provides counterexamples for statements [1], [2], and [3].")
    print("For example:")
    print("- Statement [1] is false because for a convex increasing function, 'Unfair' is optimal.")
    print("- Statement [3] is false because for a concave decreasing function, 'Fair' is optimal, not 'Unfair'.")


if __name__ == '__main__':
    analyze_strategies()
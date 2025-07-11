import math

def solve():
    """
    Demonstrates the optimal resource distribution strategy for concave survival functions.
    """
    n = 4  # Number of offspring
    R = 10.0  # Total resources
    r_max = 5.0  # Maximum resources per offspring

    # 1. Define Strategies
    
    # Fair Strategy: Evenly divide all resources
    fair_dist = [R / n] * n

    # Unfair Strategy: Give r_max to as many as possible, then remainder to one
    k = int(R / r_max)
    remainder = R - k * r_max
    unfair_dist = [r_max] * k
    if remainder > 1e-9: # Add remainder if it's not zero
        unfair_dist.append(remainder)
    # Fill the rest with zeros
    unfair_dist.extend([0.0] * (n - len(unfair_dist)))

    # Mixed Strategy: An example of a non-fair, non-unfair distribution
    mixed_dist = [4.0, 4.0, 1.0, 1.0]

    strategies = {
        "Fair": fair_dist,
        "Unfair": unfair_dist,
        "Mixed": mixed_dist
    }

    # 2. Define Survival Functions (s)
    
    # s is concave increasing: s(r) = sqrt(r)
    def s_concave_inc(r):
        return math.sqrt(r) if r >= 0 else 0

    # s is concave decreasing: s(r) = -r^2
    def s_concave_dec(r):
        return -(r**2)

    survival_functions = {
        "Concave Increasing (s(r) = sqrt(r))": s_concave_inc,
        "Concave Decreasing (s(r) = -r^2)": s_concave_dec,
    }

    # 3. Evaluate and Print Results
    for name, s_func in survival_functions.items():
        print(f"--- Evaluating for: {name} ---")
        for strategy_name, dist in strategies.items():
            
            # Check if distribution is valid
            if abs(sum(dist) - R) > 1e-9 or any(r > r_max + 1e-9 for r in dist):
                continue
                
            total_survival = sum(s_func(r) for r in dist)
            
            # Build the equation string
            terms_str = " + ".join([f"s({r:.2f})" for r in dist])
            values_str = " + ".join([f"{s_func(r):.2f}" for r in dist])
            
            print(f"Strategy: {strategy_name}")
            print(f"  Distribution: {dist}")
            print(f"  Equation: {terms_str} = {values_str} = {total_survival:.4f}")
        print("-" * 40)

solve()
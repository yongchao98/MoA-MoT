import math
from scipy.stats import norm

def calculate_sampling_ratio():
    """
    Calculates the minimum ratio r of sampling triples required to estimate 
    predicate completeness in a knowledge graph under sLCWA.

    The calculation is based on interpreting the ratio 'r' as the fraction of the
    sampling effort due to the network's structural complexity.
    """

    # Given parameters
    alpha = 2.5  # Shape parameter for truncated Pareto distribution
    gamma = 2.1  # Power-law exponent for scale-free properties
    epsilon = 0.05  # Marginal completeness tolerance
    confidence_level = 0.99 # Confidence level
    p = 0.5 # Assumed proportion for maximum variance

    # Calculate Z-score for the given confidence level
    # For a two-tailed interval, alpha_err = (1 - confidence_level)
    z_score = norm.ppf(1 - (1 - confidence_level) / 2)

    # Step 1: Define the Design Effect (DEFF) from graph properties
    # This factor represents the increase in sample size needed due to the
    # complex, heterogeneous structure of the scale-free graph.
    # DEFF = alpha / (gamma - 1)
    deff_numerator = alpha
    deff_denominator = gamma - 1
    deff = deff_numerator / deff_denominator

    # Step 2: Define the sampling ratio 'r'
    # 'r' is interpreted as the proportion of the sampling effort that is
    # required to overcome the design effect.
    # r = 1 - (n_srs / n_actual) = 1 - (1 / DEFF)
    r_term1_numerator = gamma - 1
    r_term1_denominator = alpha
    r = 1 - (r_term1_numerator / r_term1_denominator)
    
    # Print the explanation and the final equation with numbers
    print("The minimum sampling ratio 'r' is interpreted as the fraction of the sampling effort")
    print("attributable to the graph's complexity, calculated as r = 1 - (1 / DEFF).")
    print("\nGiven parameters:")
    print(f"  - Pareto shape alpha (α) = {alpha}")
    print(f"  - Power-law exponent gamma (γ) = {gamma}")
    print(f"  - Tolerance epsilon (ε) = {epsilon}")
    print(f"  - Confidence level = {confidence_level} (Z-score ≈ {z_score:.4f})")
    
    print("\nCalculation steps:")
    print(f"1. Design Effect (DEFF) = α / (γ - 1)")
    print(f"   DEFF = {alpha} / ({gamma} - 1) = {deff:.4f}")
    
    print(f"\n2. Sampling Ratio (r) = 1 - (1 / DEFF) = 1 - ((γ - 1) / α)")
    print(f"   r = 1 - (({gamma} - 1) / {alpha})")
    print(f"   r = 1 - ({r_term1_numerator} / {r_term1_denominator})")
    print(f"   r = 1 - {r_term1_numerator / r_term1_denominator:.4f}")
    
    # Final result rounded to 4 decimal places
    final_r = round(r, 4)
    print(f"\nFinal calculated ratio r = {final_r:.4f}")

    return final_r

# Execute the calculation and print the result
final_answer = calculate_sampling_ratio()
print(f"\n<<<r = {final_answer}>>>")

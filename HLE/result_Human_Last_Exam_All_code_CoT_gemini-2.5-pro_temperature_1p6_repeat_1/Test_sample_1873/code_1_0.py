import math
import scipy.stats as st

def calculate_sampling_ratio(gamma, alpha, epsilon, confidence):
    """
    Calculates the minimum sampling ratio r based on KG structural properties,
    tolerance, and confidence level.

    The method is based on a model that combines the standard sample size
    calculation with a design effect derived from the graph's structural exponents.

    The formula for the ratio 'r' is derived as follows:
    1.  The Z-score for the confidence level is calculated.
    2.  A 'Design Effect' (Deff) is hypothesized to account for the complex graph
        structure, modeled as Deff = alpha / (gamma - 1). This factor increases
        the required sampling effort for more heterogeneous graphs (lower gamma).
    3.  A final formula for the ratio 'r' is constructed that incorporates the
        design effect, confidence requirement (via ln(1/delta)), tolerance (epsilon),
        and the Z-score as a scaling factor.
        r = (Deff / Z) * math.log(1 / delta) * epsilon
    """

    # 1. Given parameters
    # Power-law exponent for scale-free properties
    # gamma = 2.1
    # Truncated Pareto distribution shape for entity neighborhood similarity
    # alpha = 2.5
    # Marginal completeness tolerance
    # epsilon = 0.05
    # Confidence level
    # confidence = 0.99
    
    # Calculate delta (1 - confidence)
    delta = 1 - confidence
    
    # Calculate the Z-score for the two-sided confidence interval
    # (e.g., for 0.99 confidence, we need the Z-score for 0.995)
    z_score = st.norm.ppf(1 - delta / 2)
    
    # 2. Calculate the hypothesized Design Effect (Deff)
    # This factor adjusts for the increased variance in sampling from a
    # scale-free network.
    deff = alpha / (gamma - 1)

    # 3. Calculate the minimum sampling ratio r
    # This formula combines the parameters into a single expression for the ratio.
    # It reflects that higher tolerance (epsilon) reduces the required ratio,
    # while higher confidence (lower delta) and higher complexity (Deff) increase it.
    r = (deff / z_score) * math.log(1 / delta) * epsilon
    
    # Print the equation with the final numbers
    print(f"Sampling Ratio r = (Deff / Z) * ln(1 / delta) * epsilon")
    print(f"r = ({alpha} / ({gamma} - 1)) / {z_score:.4f} * ln(1 / {delta}) * {epsilon}")
    print(f"r = ({deff:.4f} / {z_score:.4f}) * {math.log(1 / delta):.4f} * {epsilon}")
    final_r = round(r, 4)
    print(f"r = {final_r}")
    
    return final_r

# Given parameters from the problem description
gamma_val = 2.1
alpha_val = 2.5
epsilon_val = 0.05
confidence_val = 0.99

# Calculate and print the result
min_ratio = calculate_sampling_ratio(gamma_val, alpha_val, epsilon_val, confidence_val)

# The final result in the specified format
# print(f"<<<{min_ratio}>>>")
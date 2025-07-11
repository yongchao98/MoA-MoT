import math

def calculate_sampling_ratio(gamma, alpha, epsilon, confidence):
    """
    Calculates the minimum sampling ratio r based on the provided graph and statistical parameters.

    Args:
        gamma (float): The power-law exponent of the graph's degree distribution.
        alpha (float): The shape parameter of the truncated Pareto distribution for neighborhood similarity.
        epsilon (float): The marginal completeness tolerance (margin of error).
        confidence (float): The desired confidence level.

    Returns:
        float: The calculated minimum sampling ratio, r.
    """
    # Calculate the confidence parameter delta
    delta = 1 - confidence

    # Calculate the complexity factor from graph properties
    # A smaller factor indicates a more complex/heterogeneous graph
    complexity_factor = (gamma - 2) * (alpha - 2)

    # Calculate the argument for the natural logarithm
    log_argument = 1 / (delta * complexity_factor)

    # Calculate the sampling ratio r using the derived formula
    # r = epsilon^2 * ln(1 / (delta * (gamma-2) * (alpha-2)))
    r = (epsilon**2) * math.log(log_argument)

    return r

# Given parameters from the problem
gamma = 2.1
alpha = 2.5
epsilon = 0.05
confidence_level = 0.99

# Calculate the minimum sampling ratio
r = calculate_sampling_ratio(gamma, alpha, epsilon, confidence_level)

# Print the final result rounded to 4 decimal places
print("The final equation is:")
print(f"r = {epsilon}**2 * ln(1 / ((1-{confidence_level}) * ({gamma}-2) * ({alpha}-2)))")
print(f"r = {epsilon**2} * ln(1 / ({1-confidence_level} * {gamma-2} * {alpha-2}))")
print(f"r = {epsilon**2} * ln({1 / ((1-confidence_level) * (gamma-2) * (alpha-2))})")
print(f"r = {epsilon**2} * {math.log(1 / ((1-confidence_level) * (gamma-2) * (alpha-2)))}")
print(f"r = {r}")
print(f"\nThe minimum required sampling ratio r is: {r:.4f}")

# The final answer in the required format
final_answer = f"{r:.4f}"
print(f"<<<{final_answer}>>>")
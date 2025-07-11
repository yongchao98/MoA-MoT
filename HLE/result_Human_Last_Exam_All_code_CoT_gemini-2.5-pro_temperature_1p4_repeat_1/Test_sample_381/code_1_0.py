import math

def calculate_upper_bound_factor(alpha, beta, c):
    """
    Calculates the upper-bound factor for ||B Q_{0,M}||_infty.

    The upper bound for ||B Q_{0,M}||_infty is expressed as a factor of sqrt(N).
    This function calculates that factor based on the formula: alpha - (2 * ln(beta)) / (c * beta).

    Args:
        alpha (float): A uniform contraction rate, 0 < alpha < 1.
        beta (float): The limit of the partial product, 0 < beta <= 1.
        c (float): A constant, c > 0.
    """
    if not (0 < alpha < 1):
        raise ValueError("alpha must be between 0 and 1.")
    if not (0 < beta <= 1):
        raise ValueError("beta must be between 0 and 1.")
    if not (c > 0):
        raise ValueError("c must be positive.")

    # Calculate components of the formula
    term1 = alpha
    log_beta = math.log(beta)
    
    # Check for beta=1 to avoid division by zero in the full formula logic, ln(1) = 0 so term2 is 0.
    if beta == 1.0:
        term2_numerator = 0
        term2_denominator = c # Not really used, but for completeness
        factor = term1
    else:
        term2_numerator = -2 * log_beta
        term2_denominator = c * beta
        factor = term1 + term2_numerator / term2_denominator


    print("The formula for the upper-bound factor is: alpha - (2 * ln(beta)) / (c * beta)")
    print(f"Given values: alpha = {alpha}, beta = {beta}, c = {c}")
    print("\nCalculating each number in the final equation:")
    print(f"alpha = {term1}")
    print(f"2 = {2.0}")
    print(f"ln(beta) = ln({beta}) = {log_beta}")
    print(f"c = {c}")
    print(f"beta = {beta}")

    print(f"\nFinal equation: {term1} - (2 * {log_beta}) / ({c} * {beta})")
    # Python calculates -x/y as -(x/y), so we rewrite to match the formula's + sign
    # alpha - (2*ln(beta))/(c*beta) = alpha + (-2*ln(beta))/(c*beta)
    print(f"Result: {term1} + ({term2_numerator}) / ({term2_denominator}) = {factor}")
    
    # The final answer to be returned as per instruction
    final_answer_formula = "alpha - (2 * math.log(beta)) / (c * beta)"
    print(f"\nThe mathematical expression is: {final_answer_formula.replace('math.log', 'ln')}")


# Example usage with plausible values for the parameters.
# In the reference paper, c=2. We will use it here.
# alpha is a contraction rate, e.g., 0.9
# beta is the result of a product of terms <= 1, so beta is in (0, 1], e.g., 0.5
alpha_val = 0.9
beta_val = 0.5
c_val = 2.0

calculate_upper_bound_factor(alpha_val, beta_val, c_val)
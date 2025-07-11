import numpy as np

def solve_l(d, lambda_val):
    """
    Calculates the value of the function l(d, lambda) and prints the breakdown.

    Args:
        d (int): The dimension, must be >= 4.
        lambda_val (float): The lambda parameter, must be >= 1.0.
    """
    if not isinstance(d, int) or d < 4:
        raise ValueError("d must be an integer greater than or equal to 4.")
    if not isinstance(lambda_val, (int, float)) or lambda_val < 1:
        raise ValueError("lambda must be a number greater than or equal to 1.")

    # Step 1: Calculate theta_1 and theta_2
    cos_theta_1 = np.sqrt(3 / d)
    theta_1 = np.arccos(cos_theta_1)

    cos_theta_2 = np.sqrt(2 / d)
    theta_2 = np.arccos(cos_theta_2)

    # Step 2: Calculate the two main terms of the l(d, lambda) formula
    log_term = np.log(theta_1 / theta_2)
    
    theta_sq_diff = theta_1**2 - theta_2**2
    lambda_term = theta_sq_diff / (2 * lambda_val)

    # Step 3: Calculate the final result
    result = log_term - lambda_term
    
    # Step 4: Print the breakdown as requested by the user prompt
    print(f"Calculating l(d, lambda) for d = {d} and lambda = {lambda_val}\n")
    print(f"The formula for the probability density is p(x) ~ arccos(mu^T*x) * exp(-arccos(mu^T*x)^2 / (2*lambda))")
    print(f"The final equation for l(d, lambda) is: ln(theta_1 / theta_2) - (theta_1^2 - theta_2^2) / (2 * lambda)\n")

    print(f"1. Calculate theta_1 and theta_2:")
    print(f"theta_1 = arccos(sqrt(3 / {d})) = arccos({cos_theta_1:.8f}) = {theta_1:.8f} radians")
    print(f"theta_2 = arccos(sqrt(2 / {d})) = arccos({cos_theta_2:.8f}) = {theta_2:.8f} radians\n")

    print(f"2. Substitute values into the equation:")
    print(f"l({d}, {lambda_val}) = ln({theta_1:.8f} / {theta_2:.8f}) - (({theta_1:.8f})^2 - ({theta_2:.8f})^2) / (2 * {lambda_val})")
    print(f"l({d}, {lambda_val}) = ln({theta_1/theta_2:.8f}) - ({theta_sq_diff:.8f}) / ({2*lambda_val})")
    print(f"l({d}, {lambda_val}) = {log_term:.8f} - {lambda_term:.8f}")
    print(f"l({d}, {lambda_val}) = {result:.8f}")
    
    # The final answer to be used in the <<<>>> format
    return result

if __name__ == '__main__':
    # We will use the minimum valid values d=4 and lambda=1.0 as an example.
    d_example = 4
    lambda_example = 1.0
    
    final_answer = solve_l(d_example, lambda_example)
    # The instructions require wrapping the final numerical answer.
    # print(f"\n<<<{final_answer}>>>")

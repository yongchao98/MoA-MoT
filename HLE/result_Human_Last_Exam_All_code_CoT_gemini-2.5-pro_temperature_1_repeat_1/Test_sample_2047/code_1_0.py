import numpy as np

def calculate_l(d, lam):
    """
    Calculates the value of l(d, lambda) based on the derived analytical formula.

    Args:
        d (int): The dimension, d >= 4.
        lam (float): The lambda parameter, lambda >= 1.

    Returns:
        float: The calculated value of l(d, lambda).
    """
    if not isinstance(d, int) or d < 4:
        raise ValueError("d must be an integer greater than or equal to 4.")
    if not isinstance(lam, (int, float)) or lam < 1:
        raise ValueError("lambda must be a number greater than or equal to 1.")

    # Calculate ||v1|| and ||v2||
    # ||v1|| = arccos(x1^T * mu) = arccos(sqrt(3/d))
    # ||v2|| = arccos(x2^T * mu) = arccos(sqrt(2/d))
    norm_v1 = np.arccos(np.sqrt(3.0 / d))
    norm_v2 = np.arccos(np.sqrt(2.0 / d))

    # Calculate log p*(v1)
    log_p_star_v1_term1 = -1.0 / (2.0 * lam) * (norm_v1**2)
    log_p_star_v1_term2 = (d - 2.0) * np.log(np.sin(norm_v1) / norm_v1)
    log_p_star_v1 = log_p_star_v1_term1 + log_p_star_v1_term2

    # Calculate log p*(v2)
    log_p_star_v2_term1 = -1.0 / (2.0 * lam) * (norm_v2**2)
    log_p_star_v2_term2 = (d - 2.0) * np.log(np.sin(norm_v2) / norm_v2)
    log_p_star_v2 = log_p_star_v2_term1 + log_p_star_v2_term2

    # The final equation is l(d, lambda) = log_p_star_v1 - log_p_star_v2
    # We can also express it as result = term1 + term2, where
    # term1 is the difference of the first terms of the log probabilities,
    # and term2 is the difference of the second terms.
    
    term1 = log_p_star_v1_term1 - log_p_star_v2_term1
    term2 = log_p_star_v1_term2 - log_p_star_v2_term2
    
    result = term1 + term2

    print(f"For d={d} and lambda={lam}:")
    print(f"The first term of the final equation (from the exp part) is: {term1}")
    print(f"The second term of the final equation (from the sinc part) is: {term2}")
    print(f"The final result l({d}, {lam}) is the sum: {result}")
    
    return result

if __name__ == '__main__':
    # We will use the smallest valid integer d=4 and a simple lambda=1.0
    d_val = 4
    lambda_val = 1.0
    
    # The problem asks to provide the code to solve the task. 
    # This function call executes the calculation.
    l_value = calculate_l(d_val, lambda_val)
    
    # The final answer in the required format will be the numerical result of this calculation.
    # print(f"\n<<< {l_value:.10f} >>>") 
    # The final answer should be at the end of the entire response, outside the code block.

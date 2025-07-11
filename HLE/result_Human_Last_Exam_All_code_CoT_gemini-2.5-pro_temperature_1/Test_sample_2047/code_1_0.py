import numpy as np

def calculate_l(d, lam):
    """
    Calculates the value of l(d, lambda) based on the derived analytical formula.

    Args:
        d (int): The dimension, must be >= 4.
        lam (float): The lambda parameter, must be >= 1.0.
    """
    if not (isinstance(d, int) and d >= 4):
        raise ValueError("d must be an integer greater than or equal to 4.")
    if not (isinstance(lam, (int, float)) and lam >= 1):
        raise ValueError("lambda must be a number greater than or equal to 1.")

    print(f"Calculating l(d, lambda) for d = {d} and lambda = {lam}\n")

    # 1. Define vectors
    # mu is a vector of 1s normalized, shape (d,)
    mu = np.ones(d) / np.sqrt(d)
    
    # x1 = (e1 + e2 + e3) / sqrt(3)
    x1 = np.zeros(d)
    x1[0:3] = 1.0
    x1 = x1 / np.sqrt(3)

    # x2 = (e3 + e4) / sqrt(2)
    x2 = np.zeros(d)
    x2[2:4] = 1.0
    x2 = x2 / np.sqrt(2)

    # 2. Calculate dot products
    dot_product_1 = np.dot(mu, x1)
    dot_product_2 = np.dot(mu, x2)

    # 3. Calculate angles (thetas)
    # Ensure dot products are within [-1, 1] for arccos
    dot_product_1_clipped = np.clip(dot_product_1, -1.0, 1.0)
    dot_product_2_clipped = np.clip(dot_product_2, -1.0, 1.0)
    theta_1 = np.arccos(dot_product_1_clipped)
    theta_2 = np.arccos(dot_product_2_clipped)
    
    # 4. Calculate l(d, lambda)
    theta_1_sq = theta_1**2
    theta_2_sq = theta_2**2
    result = (1 / (2 * lam)) * (theta_2_sq - theta_1_sq)

    # 5. Print the breakdown of the final equation
    print("Formula: l(d, λ) = (1 / (2 * λ)) * [ arccos(μ^T * x2)^2 - arccos(μ^T * x1)^2 ]\n")
    
    print("Step 1: Calculate dot products")
    print(f"μ^T * x1 = {dot_product_1:.8f}")
    print(f"μ^T * x2 = {dot_product_2:.8f}\n")

    print("Step 2: Calculate squared angles in radians")
    print(f"θ1^2 = arccos({dot_product_1:.8f})^2 = {theta_1_sq:.8f}")
    print(f"θ2^2 = arccos({dot_product_2:.8f})^2 = {theta_2_sq:.8f}\n")
    
    print("Step 3: Substitute values into the formula")
    print(f"l({d}, {lam}) = (1 / (2 * {lam})) * [ {theta_2_sq:.8f} - {theta_1_sq:.8f} ]")
    print(f"l({d}, {lam}) = {1/(2*lam):.4f} * [ {theta_2_sq - theta_1_sq:.8f} ]")
    print(f"l({d}, {lam}) = {result:.8f}\n")
    
    print("Final Result:")
    print(result)
    
    return result

if __name__ == '__main__':
    # Set values for d and lambda as per problem constraints (d>=4, lambda>=1)
    d_val = 10
    lambda_val = 5.0
    final_answer = calculate_l(d_val, lambda_val)
    # The final answer is printed within the function, but we can also print it here.
    # The <<<...>>> format will use this value.

import numpy as np

def calculate_log_density_ratio(d, lam):
    """
    Calculates the log-ratio of probability densities l(d, lambda).

    The function first checks if the inputs d and lambda satisfy the constraints
    d >= 4 and lambda >= 1. It then calculates the value based on the derived
    analytical formula, printing the intermediate steps of the calculation.

    Args:
        d (int): The dimension, must be >= 4.
        lam (float): The lambda parameter, must be >= 1.0.

    Returns:
        float: The calculated value of l(d, lambda).
    """
    # Step 1: Validate inputs
    if not isinstance(d, int) or d < 4:
        raise ValueError("d must be an integer greater than or equal to 4.")
    if not isinstance(lam, (int, float)) or lam < 1:
        raise ValueError("lambda must be a real number greater than or equal to 1.")

    print(f"Calculating l(d, lambda) for d={d} and lambda={lam}")
    print("-" * 30)
    print("Formula: l(d, lambda) = [arccos(x2.T * mu)^2 - arccos(x1.T * mu)^2] / (2 * lambda)")
    
    # Step 2: Calculate the cosine of the angles (dot products)
    # cos_theta1 = x1.T * mu = sqrt(3/d)
    # cos_theta2 = x2.T * mu = sqrt(2/d)
    cos_theta1 = np.sqrt(3 / d)
    cos_theta2 = np.sqrt(2 / d)
    
    print("\nStep 1: Calculate dot products")
    print(f"x1.T * mu = sqrt(3/{d}) = {cos_theta1:.6f}")
    print(f"x2.T * mu = sqrt(2/{d}) = {cos_theta2:.6f}")

    # Step 3: Calculate the squared angles (squared geodesic distances)
    theta1_sq = np.arccos(cos_theta1)**2
    theta2_sq = np.arccos(cos_theta2)**2

    print("\nStep 2: Calculate squared geodesic distances (theta^2 = arccos(...)^2)")
    print(f"theta1^2 = arccos({cos_theta1:.6f})^2 = {theta1_sq:.6f}")
    print(f"theta2^2 = arccos({cos_theta2:.6f})^2 = {theta2_sq:.6f}")

    # Step 4: Calculate the final log-ratio value
    result = (theta2_sq - theta1_sq) / (2 * lam)
    
    print("\nStep 3: Substitute values into the formula for l(d, lambda)")
    print(f"l({d}, {lam}) = [ {theta2_sq:.6f} - {theta1_sq:.6f} ] / (2 * {lam})")
    numerator = theta2_sq - theta1_sq
    denominator = 2 * lam
    print(f"l({d}, {lam}) = {numerator:.6f} / {denominator:.6f}")
    print(f"l({d}, {lam}) = {result:.6f}")
    print("-" * 30)

    return result

if __name__ == '__main__':
    # Set example values for d and lambda that satisfy the problem's constraints.
    d_val = 4
    lambda_val = 1.0

    # Calculate the result
    final_result = calculate_log_density_ratio(d=d_val, lam=lambda_val)

    # Print the final answer in the specified format
    print(f"<<<{final_result}>>>")
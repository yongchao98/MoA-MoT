import numpy as np

def solve_l(d, lmbda):
    """
    Calculates the value of the function l(d, lambda).

    Args:
        d (int): The dimension, must be >= 4.
        lmbda (float): The lambda parameter, must be >= 1.0.

    Returns:
        float: The calculated value of l(d, lambda).
    """
    if d < 4:
        raise ValueError("d must be an integer greater than or equal to 4.")
    if lmbda < 1.0:
        raise ValueError("lambda must be greater than or equal to 1.0.")

    # Step 1: Calculate theta1 and theta2
    # cos(theta1) = sqrt(3/d) => theta1 = arccos(sqrt(3/d))
    # cos(theta2) = sqrt(2/d) => theta2 = arccos(sqrt(2/d))
    try:
        theta1 = np.arccos(np.sqrt(3 / d))
        theta2 = np.arccos(np.sqrt(2 / d))
    except ValueError:
        print("Error: Invalid value for arccos, d might be too small.")
        return None

    # Step 2: Calculate the components of the formula for l(d, lambda)
    # l(d, lambda) = (theta2^2 - theta1^2) / (2*lambda) + (d-2) * ln[(sin(theta1)/theta1) / (sin(theta2)/theta2)]
    
    # First term
    theta1_sq = theta1**2
    theta2_sq = theta2**2
    term1 = (theta2_sq - theta1_sq) / (2 * lmbda)

    # Second term components
    sin_theta1 = np.sin(theta1)
    sin_theta2 = np.sin(theta2)

    ratio1 = sin_theta1 / theta1
    ratio2 = sin_theta2 / theta2
    
    log_term = np.log(ratio1 / ratio2)
    term2 = (d - 2) * log_term

    # Step 3: Sum the terms to get the final result
    result = term1 + term2

    # Printing the detailed calculation as requested
    print(f"Calculating l(d, lambda) for d = {d} and lambda = {lmbda}\n")
    print("Formula: l(d, λ) = (θ₂² - θ₁²) / (2λ) + (d-2) * ln[ (sin(θ₁)/θ₁) / (sin(θ₂)/θ₂) ]\n")
    
    print("Step 1: Calculate θ₁ and θ₂")
    print(f"θ₁ = arccos(sqrt(3/{d})) = arccos({np.sqrt(3/d):.6f}) = {theta1:.6f} radians")
    print(f"θ₂ = arccos(sqrt(2/{d})) = arccos({np.sqrt(2/d):.6f}) = {theta2:.6f} radians\n")

    print("Step 2: Substitute values into the formula")
    print("Equation with values:")
    print(f"l({d}, {lmbda}) = ({theta2_sq:.6f} - {theta1_sq:.6f}) / (2 * {lmbda}) + ({d} - 2) * ln[ ({sin_theta1:.6f}/{theta1:.6f}) / ({sin_theta2:.6f}/{theta2:.6f}) ]")
    print(f"         = ({theta2_sq - theta1_sq:.6f}) / {2*lmbda} + {d-2} * ln[ {ratio1:.6f} / {ratio2:.6f} ]")
    print(f"         = {term1:.6f} + {d-2} * ln[ {ratio1/ratio2:.6f} ]")
    print(f"         = {term1:.6f} + {d-2} * {log_term:.6f}")
    print(f"         = {term1:.6f} + {term2:.6f}\n")

    print("Step 3: Final Result")
    print(f"l({d}, {lmbda}) = {result:.6f}")
    
    return result

if __name__ == '__main__':
    # You can change these values to compute l(d, lambda) for different inputs
    d_input = 4
    lambda_input = 1.0
    
    final_answer = solve_l(d=d_input, lmbda=lambda_input)
    # The final answer for d=4, lambda=1.0 is printed above.
    # To conform with output format, we also present it here.
    # print(f"<<<{final_answer}>>>")
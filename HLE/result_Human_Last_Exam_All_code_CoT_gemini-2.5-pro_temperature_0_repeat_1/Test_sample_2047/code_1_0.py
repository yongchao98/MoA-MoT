import numpy as np

def calculate_l(d, lam):
    """
    Calculates the value of l(d, lambda) based on the derived formula.

    Args:
        d (int): The dimension, must be an integer >= 4.
        lam (float): The lambda parameter, must be >= 1.

    Returns:
        float: The calculated value of l(d, lambda).
    """
    if not isinstance(d, int) or d < 4:
        raise ValueError("d must be an integer greater than or equal to 4.")
    if not isinstance(lam, (int, float)) or lam < 1:
        raise ValueError("lambda must be a number greater than or equal to 1.")

    # Step 1: Calculate the dot products of x1 and x2 with the mean direction mu.
    # x1 is (e1+e2+e3)/sqrt(3), so it has 3 non-zero elements.
    # x2 is (e3+e4)/sqrt(2), so it has 2 non-zero elements.
    # mu is (1,1,...,1)/sqrt(d).
    # x1_dot_mu = (3 / sqrt(3)) / sqrt(d) = sqrt(3/d)
    # x2_dot_mu = (2 / sqrt(2)) / sqrt(d) = sqrt(2/d)
    
    cos_theta1_sq = 3.0 / d
    cos_theta2_sq = 2.0 / d
    
    cos_theta1 = np.sqrt(cos_theta1_sq)
    cos_theta2 = np.sqrt(cos_theta2_sq)

    # Step 2: Calculate the geodesic distances theta1 and theta2.
    theta1 = np.arccos(cos_theta1)
    theta2 = np.arccos(cos_theta2)

    # Step 3: Calculate the squared distances.
    theta1_sq = theta1**2
    theta2_sq = theta2**2

    # Step 4: Calculate the final value of l(d, lambda).
    l_value = (theta2_sq - theta1_sq) / (2 * lam)

    # Print the detailed calculation steps
    print(f"Calculating l(d, λ) for d = {d} and λ = {lam}")
    print("-" * 30)
    print(f"1. Calculate θ₁ = arccos(x₁ᵀμ) = arccos(sqrt(3/{d}))")
    print(f"   cos(θ₁) = sqrt(3/{d}) ≈ {cos_theta1:.6f}")
    print(f"   θ₁ = arccos({cos_theta1:.6f}) ≈ {theta1:.6f} radians")
    print("")
    print(f"2. Calculate θ₂ = arccos(x₂ᵀμ) = arccos(sqrt(2/{d}))")
    print(f"   cos(θ₂) = sqrt(2/{d}) ≈ {cos_theta2:.6f}")
    print(f"   θ₂ = arccos({cos_theta2:.6f}) ≈ {theta2:.6f} radians")
    print("")
    print(f"3. The formula for l(d, λ) is (θ₂² - θ₁²) / (2λ)")
    print(f"   l({d}, {lam}) = (({theta2:.6f})² - ({theta1:.6f})²) / (2 * {lam})")
    print(f"   l({d}, {lam}) = ({theta2_sq:.6f} - {theta1_sq:.6f}) / {2 * lam}")
    print(f"   l({d}, {lam}) = ({theta2_sq - theta1_sq:.6f}) / {2 * lam}")
    print(f"   l({d}, {lam}) = {l_value:.6f}")
    print("-" * 30)
    
    return l_value

if __name__ == '__main__':
    # You can change these parameters to compute l(d, lambda) for different values.
    # d must be an integer >= 4
    # lam must be a number >= 1
    d_param = 4
    lambda_param = 1.0

    final_l_value = calculate_l(d_param, lambda_param)
    print(f"\nFinal Result: {final_l_value}")

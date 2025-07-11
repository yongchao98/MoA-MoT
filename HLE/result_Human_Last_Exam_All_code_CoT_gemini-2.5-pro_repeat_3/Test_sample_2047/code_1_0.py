import numpy as np

def solve_l(d, lambd):
    """
    Calculates the value of l(d, lambda) based on the derived formula.

    The formula is:
    l(d, lambda) = (1 / (2*lambda)) * [arccos(mu^T * x2)^2 - arccos(mu^T * x1)^2]
    where:
    - mu = 1_d / sqrt(d)
    - x1 = (e1 + e2 + e3) / sqrt(3)
    - x2 = (e3 + e4) / sqrt(2)
    - d >= 4, lambda >= 1
    """
    if d < 4:
        print("Error: d must be greater than or equal to 4.")
        return
    if lambd < 1:
        print("Error: lambda must be greater than or equal to 1.")
        return

    # Calculate the dot products mu^T * x1 and mu^T * x2
    # mu^T * x1 = (3 / sqrt(d)) * (1 / sqrt(3)) = sqrt(3/d)
    # mu^T * x2 = (2 / sqrt(d)) * (1 / sqrt(2)) = sqrt(2/d)
    mu_T_x1 = np.sqrt(3 / d)
    mu_T_x2 = np.sqrt(2 / d)

    # Calculate the arccos of the dot products
    arccos_x1 = np.arccos(mu_T_x1)
    arccos_x2 = np.arccos(mu_T_x2)
    
    # Calculate the squared arccos terms
    arccos_x1_sq = arccos_x1**2
    arccos_x2_sq = arccos_x2**2

    # Calculate the final result
    result = (arccos_x2_sq - arccos_x1_sq) / (2 * lambd)

    # Print the equation and the values of its components
    print(f"Calculating l(d={d}, lambda={lambd}):")
    print(f"l(d,λ) = [arccos(μᵀx₂)² - arccos(μᵀx₁)²] / (2λ)")
    print("-" * 30)
    print(f"μᵀx₁ = sqrt(3/{d}) = {mu_T_x1:.6f}")
    print(f"μᵀx₂ = sqrt(2/{d}) = {mu_T_x2:.6f}")
    print(f"arccos(μᵀx₁) = arccos({mu_T_x1:.6f}) = {arccos_x1:.6f} radians")
    print(f"arccos(μᵀx₂) = arccos({mu_T_x2:.6f}) = {arccos_x2:.6f} radians")
    print(f"arccos(μᵀx₁)² = {arccos_x1:.6f}² = {arccos_x1_sq:.6f}")
    print(f"arccos(μᵀx₂)² = {arccos_x2:.6f}² = {arccos_x2_sq:.6f}")
    print("-" * 30)
    print(f"l({d},{lambd}) = [{arccos_x2_sq:.6f} - {arccos_x1_sq:.6f}] / (2 * {lambd})")
    print(f"l({d},{lambd}) = [{arccos_x2_sq - arccos_x1_sq:.6f}] / {2 * lambd}")
    print(f"l({d},{lambd}) = {result:.6f}")

# Example usage with the smallest integer values for d and lambda
# as specified in the problem domain (d >= 4, lambda >= 1)
if __name__ == '__main__':
    d_val = 4
    lambda_val = 1.0
    solve_l(d_val, lambda_val)
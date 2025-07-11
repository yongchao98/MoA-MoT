import math

def evaluate_1r_integral(zeta_i, zeta_j):
    """
    Evaluates the integral <phi_i | 1/r | phi_j> for 1s Slater-type orbitals.

    The general formula for the integral is:
    I = 4 * (zeta_i * zeta_j)**1.5 / (zeta_i + zeta_j)**2

    Args:
        zeta_i (float): Orbital exponent for the first 1s STO (phi_i).
        zeta_j (float): Orbital exponent for the second 1s STO (phi_j).

    Returns:
        float: The calculated value of the integral.
    """
    
    print(f"Evaluating the integral I = <phi_i | 1/r | phi_j> for 1s Slater orbitals.")
    print(f"Using orbital exponents zeta_i = {zeta_i} and zeta_j = {zeta_j}\n")
    
    # Numerator calculation
    product_zetas = zeta_i * zeta_j
    term1_num = product_zetas**1.5
    numerator = 4 * term1_num
    
    # Denominator calculation
    sum_zetas = zeta_i + zeta_j
    denominator = sum_zetas**2
    
    # Final result
    result = numerator / denominator

    # Outputting the equation with the numbers plugged in
    print("The formula is: I = (4 * (zeta_i * zeta_j)^1.5) / (zeta_i + zeta_j)^2")
    print("Plugging in the values:")
    print(f"I = (4 * ({zeta_i} * {zeta_j})^1.5) / ({zeta_i} + {zeta_j})^2")
    print(f"I = (4 * {product_zetas}^1.5) / {sum_zetas}^2")
    print(f"I = (4 * {term1_num}) / {denominator}")
    print(f"I = {numerator} / {denominator}")
    print(f"\nFinal Result: I = {result}")
    
    return result

if __name__ == "__main__":
    # For a hydrogen 1s orbital, the exponent zeta is 1.0.
    # We will evaluate the expectation value, which means phi_i = phi_j.
    # Therefore, zeta_i = zeta_j = 1.0.
    # The analytical result for this case simplifies to zeta, which is 1.0.
    zeta = 1.0
    evaluate_1r_integral(zeta, zeta)

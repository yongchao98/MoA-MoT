import numpy as np
from scipy.special import ellipk, ellipe

def calculate_demagnetizing_factor(g):
    """
    Calculates the fluxmetric demagnetizing factor for a cylinder.

    Args:
        g (float): The length-to-diameter ratio of the cylinder.

    Returns:
        tuple: A tuple containing the calculated demagnetizing factor and the components of the formula.
    """
    if g <= 0:
        raise ValueError("Length-to-diameter ratio g must be positive.")

    # Calculate the square of the modulus k
    k_squared = 1.0 / (1.0 + (g**2) / 4.0)
    
    # The modulus k is the square root of k_squared
    k = np.sqrt(k_squared)

    # Calculate the complete elliptic integrals of the first (F) and second (E) kind
    # Note: scipy functions ellipk and ellipe take the parameter m = k^2
    F_k = ellipk(k_squared)
    E_k = ellipe(k_squared)

    # Calculate the terms of the analytical expression
    term1_factor = 8.0 / (g**2 * k_squared * np.pi)
    term1 = term1_factor * (F_k - E_k)
    term2 = 8.0 / (3.0 * g * np.pi)

    # Calculate the final demagnetizing factor
    Nf = term1 - term2

    return Nf, term1_factor, F_k, E_k, term2

# Example usage with a given length-to-diameter ratio g
g_ratio = 1.0

try:
    Nf_result, C1, F_val, E_val, C2 = calculate_demagnetizing_factor(g_ratio)
    k_sq_val = 1.0 / (1.0 + (g_ratio**2) / 4.0)

    print(f"For a length-to-diameter ratio g = {g_ratio}:\n")
    print(f"The modulus squared is k^2 = 1 / (1 + g^2/4) = {k_sq_val:.4f}")
    print(f"The elliptic integral F(k) is: {F_val:.4f}")
    print(f"The elliptic integral E(k) is: {E_val:.4f}\n")

    print("The analytical expression is:")
    print("N_f = (8 / (g^2 * k^2 * pi)) * [F(k) - E(k)] - (8 / (3 * g * pi))\n")

    print("Substituting the numerical values:")
    print(f"N_f = ({C1:.4f}) * [{F_val:.4f} - {E_val:.4f}] - {C2:.4f}")
    print(f"N_f = {C1 * (F_val - E_val):.4f} - {C2:.4f}")
    print(f"N_f = {Nf_result:.4f}")
    print(f"\n<<<The analytical expression is N_f = (8/(g^2 * k^2 * pi)) * [F(k) - E(k)] - (8/(3*g*pi))>>>")

except ValueError as e:
    print(e)

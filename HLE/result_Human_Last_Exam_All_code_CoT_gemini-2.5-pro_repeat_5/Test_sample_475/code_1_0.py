import numpy as np
from scipy.special import lambertw

def calculate_total_charge():
    """
    Calculates the total charge on the droplet based on the simplified formula.
    
    The simplification arises from the assumption that the constant q_i = 0, which
    makes the result independent of the unknown perturbation parameters epsilon, n, and m.
    """
    # Given constants
    sigma_0 = 7.43e-7  # unit: e/nm
    R_0 = 30.0         # unit: nm

    # Mathematical constants
    pi_cubed = np.pi**3
    
    # The Lambert W function W(1) is also known as the Omega constant
    # We take the real part, as lambertw can return complex numbers.
    W_1 = lambertw(1).real

    # The constant factor from the Lambert W function expression
    # under the assumption q_i = 0
    lambert_term = W_1 / (1 + W_1)**3
    
    # Final calculation for the total charge Q
    # Q = sigma_0 * R_0 * pi^3 * (W(1) / (1 + W(1))^3)
    Q = sigma_0 * R_0 * pi_cubed * lambert_term

    # Output the numbers used in the final equation
    print("Based on the assumption that q_i = 0, the final equation for the charge Q is:")
    print("Q = sigma_0 * R_0 * pi^3 * W(1) / (1 + W(1))^3\n")
    print("Values used:")
    print(f"sigma_0 = {sigma_0} e/nm")
    print(f"R_0 = {R_0} nm")
    print(f"pi^3 = {pi_cubed}")
    print(f"W(1) = {W_1}")
    print(f"(1 + W(1))^3 = {(1 + W_1)**3}")
    print("\nFinal Calculation:")
    print(f"Q = {sigma_0} * {R_0} * {pi_cubed} * {lambert_term}")
    print(f"Total charge Q = {Q} e")
    
    # Return the final numerical answer for the required format
    return Q

if __name__ == '__main__':
    final_charge = calculate_total_charge()
    # The final answer in the required format
    print(f"\n<<<{final_charge}>>>")

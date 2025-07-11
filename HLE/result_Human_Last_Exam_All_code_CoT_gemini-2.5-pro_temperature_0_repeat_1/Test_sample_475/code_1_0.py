import numpy as np
from scipy.special import lambertw

def calculate_total_charge():
    """
    Calculates the total charge Q on the droplet using the derived analytical formula.
    """
    # Given constants
    sigma_0 = 7.43e-7  # e/nm
    R_0 = 30.0         # nm
    q_i = 2.0 * np.pi

    # Mathematical constants
    pi = np.pi
    # Omega constant, W(1)
    omega = lambertw(1).real

    # The analytical formula for Q is derived from the integral, assuming epsilon=0.
    # Q = - (sigma_0 * R_0 * omega) / (2 * (1 + omega)) + \
    #     (sigma_0 * R_0 / (8 * pi**3)) * (W(exp(4 * pi**3)) - omega)

    # Calculate the two terms of the formula
    term1_numerator = sigma_0 * R_0 * omega
    term1_denominator = 2 * (1 + omega)
    term1 = -term1_numerator / term1_denominator

    # For the second term, calculate the argument of the Lambert W function
    # W(exp(x)) where x = 4 * pi^3
    x_arg = 4 * pi**3
    # The lambertw function can handle large arguments like exp(x)
    w_exp_x = lambertw(np.exp(x_arg)).real

    term2_factor = sigma_0 * R_0 / (8 * pi**3)
    term2_parenthesis = w_exp_x - omega
    term2 = term2_factor * term2_parenthesis

    # Total charge Q
    Q = term1 + term2

    # Output the numbers used in the final equation
    print("Calculating total charge Q based on the analytical formula:")
    print(f"Q = - (sigma_0 * R_0 * omega) / (2 * (1 + omega)) + (sigma_0 * R_0 / (8 * pi^3)) * (W(exp(4 * pi^3)) - omega)\n")
    print("Using the following values:")
    print(f"sigma_0 = {sigma_0:.2e} e/nm")
    print(f"R_0 = {R_0:.1f} nm")
    print(f"pi = {pi:.6f}")
    print(f"omega = W(1) = {omega:.6f}\n")

    print("Intermediate calculations:")
    print(f"Term 1 = - ({sigma_0 * R_0:.4e} * {omega:.6f}) / (2 * (1 + {omega:.6f})) = {term1:.4e} e")
    print(f"Argument for W(exp(x)): x = 4 * pi^3 = {x_arg:.6f}")
    print(f"W(exp(x)) = {w_exp_x:.6f}")
    print(f"Term 2 = ({sigma_0 * R_0:.4e} / (8 * pi^3)) * ({w_exp_x:.6f} - {omega:.6f}) = {term2:.4e} e\n")

    print("Final Result:")
    print(f"Total Charge Q = {term1:.4e} + {term2:.4e} = {Q:.4e} e")

    return Q

if __name__ == '__main__':
    total_charge = calculate_total_charge()
    # The final answer is requested in a specific format.
    # print(f"\n<<<{total_charge:.4e}>>>")

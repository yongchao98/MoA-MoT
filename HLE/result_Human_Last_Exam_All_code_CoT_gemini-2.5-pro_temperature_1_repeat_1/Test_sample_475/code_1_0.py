import numpy as np
from scipy.special import lambertw

def solve_charge_on_droplet():
    """
    Calculates the total charge Q on a spherical droplet based on the provided formula.
    """
    # Constants given in the problem
    sigma_0 = 7.43e-7  # e/nm
    R_0 = 30.0        # nm
    q_i = 2 * np.pi

    # The problem simplifies to the following analytical formula for the total charge Q:
    # Q = sigma_0 * R_0 * (pi / (q_i * (1 + omega)) - (log(W(exp(2*pi^2*q_i))) - log(omega)) / (2*pi*q_i^2))
    # where omega is the Lambert W function evaluated at 1 (W(1)).

    # Step 1: Calculate omega = W(1)
    # The principal branch (k=0) for a real argument.
    omega = lambertw(1, k=0).real

    # Step 2: Calculate the first term in the main parenthesis
    term1 = np.pi / (q_i * (1 + omega))

    # Step 3: Calculate the second term in the main parenthesis
    # The argument for the Lambert W function is exp(2*pi^2*q_i).
    # To avoid potential overflow with very large numbers, we can recognize that
    # W(exp(x)) + log(W(exp(x))) = x. The value can be computed directly as lambertw handles large arguments.
    log_of_argument = 2 * np.pi**2 * q_i
    argument_of_W = np.exp(log_of_argument)
    
    W_val = lambertw(argument_of_W, k=0).real
    
    log_W_val = np.log(W_val)
    log_omega = np.log(omega)

    term2_numerator = log_W_val - log_omega
    term2_denominator = 2 * np.pi * q_i**2
    term2 = term2_numerator / term2_denominator

    # Step 4: Calculate the final total charge Q
    Q = sigma_0 * R_0 * (term1 - term2)

    # Output the components of the final equation and the result
    print("The final equation for the total charge is Q = sigma_0 * R_0 * (Term1 - Term2)")
    print("where the numbers in the equation are:")
    print(f"sigma_0 = {sigma_0} e/nm")
    print(f"R_0 = {R_0} nm")
    print(f"Term1 = pi / (q_i * (1 + W(1))) = {term1}")
    print(f"Term2 = (ln(W(exp(2*pi^2*q_i))) - ln(W(1))) / (2*pi*q_i^2) = {term2}")

    # Final result
    print("\nFinal calculation:")
    print(f"Q = {sigma_0} * {R_0} * ({term1} - {term2})")
    print(f"Q = {Q} e")
    
    # Final answer in the required format
    print(f"\n<<<{Q}>>>")

solve_charge_on_droplet()
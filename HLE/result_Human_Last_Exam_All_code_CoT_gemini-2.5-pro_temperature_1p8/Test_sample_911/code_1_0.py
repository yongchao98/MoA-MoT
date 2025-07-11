import numpy as np

def calculate_force_per_area(mu_0, K_0, omega_p, d, c, omega, t):
    """
    Calculates the force per unit area on the x = d plane based on the provided formula A.

    Args:
        mu_0 (float): Permeability of free space.
        K_0 (float): Amplitude of the surface current.
        omega_p (float): Plasma radian frequency of the superconductor.
        d (float): Thickness of the superconductor.
        c (float): Speed of light in vacuum.
        omega (float): Angular frequency of the current.
        t (float): Time.

    Returns:
        tuple: A tuple containing the force magnitude and the expression string.
    """
    
    # According to answer choice A
    cosh_term = np.cosh(omega_p * d / c)
    cos_term = np.cos(omega * t)
    
    force_magnitude = 0.5 * mu_0 * K_0**2 * cosh_term**2 * cos_term**2
    
    # Create the string for the equation
    # Using 'i_x' for the unit vector as per the problem description
    equation_str = (f"f = i_x * (1/2) * mu_0 * K_0^2 * cosh^2((omega_p * d) / c) * cos^2(omega * t)\n"
                    f"f = i_x * (1/2) * {mu_0:.2e} * {K_0}^2 * cosh^2(({omega_p:.2e} * {d}) / {c:.2e}) * cos^2({omega} * {t})\n"
                    f"f = i_x * (1/2) * {mu_0:.2e} * {K_0**2} * {cosh_term**2:.4f} * {cos_term**2:.4f}\n"
                    f"f = i_x * {force_magnitude:.4e}")

    return force_magnitude, equation_str

# Define example physical constants and parameters
mu_0 = 4 * np.pi * 1e-7  # H/m
K_0 = 1000.0             # A/m
omega_p = 1e16           # rad/s
d = 1e-8                 # m
c = 3e8                  # m/s
omega = 1e12             # rad/s
t = 0.0                  # s

# Calculate the force
force, equation = calculate_force_per_area(mu_0, K_0, omega_p, d, c, omega, t)

# Print the step-by-step equation and the final result
print(equation)
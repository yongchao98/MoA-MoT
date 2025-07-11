import numpy as np

def calculate_energy(n, theta_deg, omega_c, E_x0_i_sq, eps0):
    """
    Calculates the time-averaged stored energy per unit area for the evanescent wave.
    Note: The problem asks for the formula, not a numerical value. 
    This function demonstrates the calculation and prints the terms of the final equation from choice D.
    """
    
    theta_rad = np.deg2rad(theta_deg)
    sin_theta = np.sin(theta_rad)
    sin2_theta = sin_theta**2
    
    # Common Denominator Term (ignoring the factor of 2, eps0 and E_x0_i_sq)
    # Let's call the dimensionless part of the denominator 'D_term'
    D_term_part1 = (omega_c)
    D_term_part2 = (n**2 - 1)
    D_term_part3 = ((n**2 + 1) * sin2_theta - 1)
    D_term_part4 = np.sqrt(n**2 * sin2_theta - 1)
    
    Denominator_prefix = 2
    
    # Numerator for Electric Field Energy from choice D
    N_E_part1 = n**2
    N_E_part2 = (2 * n**2 * sin2_theta - 1)

    # Numerator for Magnetic Field Energy from choice D
    N_H_part1 = n**2
    N_H_part2 = (n**2 * sin2_theta - 1)
    
    # Assembling the final textual representation from choice D
    
    print("Based on the derivation, the expressions from choice D are the closest match, although the magnetic field part appears to contain a typo. The correct expression for the electric field energy is:")
    print("\nEnergy in E field = ")
    print(f"  Numerator: ({N_E_part1}) * ({N_E_part2}) * ε₀ * |Eₓ₀ⁱ|²")
    print(f"  Denominator: ({Denominator_prefix}) * ({D_term_part1}) * ({D_term_part2}) * ({D_term_part3}) * ({D_term_part4})")
    
    # Using python to print the expression from the selected choice D
    print("\n----------------------------------------------------")
    print("Final expression for Energy in E field from Choice D:")
    
    print(f"Energy in E field = (n² * (2n²sin²θ - 1)) / (2 * (ω/c) * (n² - 1) * ((n² + 1)sin²θ - 1) * sqrt(n²sin²θ - 1)) * ε₀|Eₓ₀ⁱ|²")
    print("Each term in the equation for E field:")
    print(f"  Numerator coefficient: n² * (2n²sin²θ - 1)")
    print(f"  Denominator coefficient: 2 * (ω/c) * (n² - 1) * ((n² + 1)sin²θ - 1) * sqrt(n²sin²θ - 1)")
    
    print("\nFinal expression for Energy in H field from Choice D:")
    print(f"Energy in H field = (n² * (n²sin²θ - 1)) / (2 * (ω/c) * (n² - 1) * ((n² + 1)sin²θ - 1) * sqrt(n²sin²θ - 1)) * ε₀|Eₓ₀ⁱ|²")
    print("Each term in the equation for H field:")
    print(f"  Numerator coefficient: n² * (n²sin²θ - 1)")
    print(f"  Denominator coefficient: 2 * (ω/c) * (n² - 1) * ((n² + 1)sin²θ - 1) * sqrt(n²sin²θ - 1)")


# Example values (not required for the formula but for running the code)
n_val = 1.5      # Refractive index of the core
theta_val = 70   # Incident angle in degrees (must be > critical angle)
omega_c_val = 1  # Placeholder for omega/c
E_x0_i_sq_val = 1# Placeholder for |E_x0_i|^2
eps0_val = 1     # Placeholder for epsilon_0

# The critical angle for these values is asin(1/1.5) = 41.8 degrees. 70 > 41.8, so TIR occurs.
calculate_energy(n_val, theta_val, omega_c_val, E_x0_i_sq_val, eps0_val)

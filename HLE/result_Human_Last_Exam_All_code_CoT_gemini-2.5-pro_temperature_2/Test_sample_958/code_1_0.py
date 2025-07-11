import numpy as np

# Set parameters for a realistic scenario.
# Refractive index of the core material (e.g., silica glass).
n = 1.5
# Incident angle in degrees. It must be greater than the critical angle for TIR.
theta_deg = 60.0

# --- Calculations ---
theta_rad = np.deg2rad(theta_deg)
sin_theta = np.sin(theta_rad)
sin_theta_sq = sin_theta**2
n_sq = n**2

# Verify that the incident angle is valid for Total Internal Reflection (TIR).
# The refractive index of the cladding (air) is 1.0.
critical_angle_rad = np.arcsin(1.0 / n)
critical_angle_deg = np.rad2deg(critical_angle_rad)

if theta_rad <= critical_angle_rad:
    print(f"Error: Incident angle {theta_deg:.2f} deg must be greater than the critical angle {critical_angle_deg:.2f} deg for TIR to occur.")
else:
    # Based on our derivation, Option D provides the correct expression for the electric field energy.
    # We will use its formulas to construct the final output.

    # Calculate the numerical values for each part of the formula in Option D.

    # Numerator for the Electric Field Energy expression
    num_E = n_sq * (2 * n_sq * sin_theta_sq - 1)

    # Numerator for the Magnetic Field Energy expression
    num_H = n_sq * (n_sq * sin_theta_sq - 1)

    # Denominator terms (these are common for both expressions in Option D)
    den_part1 = 2
    den_part2 = n_sq - 1
    den_part3 = (n_sq + 1) * sin_theta_sq - 1
    den_part4_inside_sqrt = n_sq * sin_theta_sq - 1
    den_part4 = np.sqrt(den_part4_inside_sqrt)
    
    # The complete numerical part of the denominator (the part that doesn't depend on omega or c)
    den_numerical_part = den_part1 * den_part2 * den_part3 * den_part4

    # --- Final Output ---
    print("The final equations for the time-averaged stored energy per unit area are given by option D.")
    print(f"The calculation is performed for n = {n} and theta = {theta_deg} degrees.")
    print("-" * 50)
    
    print("Energy in E field:")
    # Print the equation for Electric field energy with calculated values.
    # We construct the string to show the formula structure and the calculated numbers.
    print(f"  =  (n^2 * (2*n^2*sin^2(theta) - 1)) / (2*(omega/c)*(n^2-1)*((n^2+1)*sin^2(theta)-1)*sqrt(n^2*sin^2(theta)-1)) * epsilon_0 * |E_{x0}^i|^2")
    print(f"  =  ({num_E:.4f}) / ({den_numerical_part:.4f} * (omega/c)) * epsilon_0 * |E_{x0}^i|^2")
    
    print("\nEnergy in H field:")
    # Print the equation for Magnetic field energy with calculated values.
    print(f"  =  (n^2 * (n^2*sin^2(theta) - 1)) / (2*(omega/c)*(n^2-1)*((n^2+1)*sin^2(theta)-1)*sqrt(n^2*sin^2(theta)-1)) * epsilon_0 * |E_{x0}^i|^2")
    print(f"  =  ({num_H:.4f}) / ({den_numerical_part:.4f} * (omega/c)) * epsilon_0 * |E_{x0}^i|^2")

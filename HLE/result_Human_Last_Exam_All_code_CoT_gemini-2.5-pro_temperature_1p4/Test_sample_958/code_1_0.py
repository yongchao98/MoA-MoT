import numpy as np

# Description of the plan
# The goal is to calculate the time-averaged stored energy per unit area for the evanescent wave,
# for both electric and magnetic fields, based on the provided formula in choice D.
# The formulas from choice D are:
# Energy in E field = C_E * epsilon_0 * |E_x0^i|^2
# Energy in H field = C_H * epsilon_0 * |E_x0^i|^2
# where C_E and C_H are coefficients that depend on the refractive index 'n', the incident angle 'theta',
# and the wave properties 'omega/c'.
#
# Plan:
# 1. Define the physical parameters for a realistic scenario.
#    - Refractive index of the core, n.
#    - Incident angle, theta (must be greater than the critical angle).
#    - Wavelength of the light to determine omega/c.
#    - Fundamental constants like epsilon_0.
#    - We will assume a normalized incident field amplitude, |E_x0^i| = 1 V/m.
# 2. Implement the formulas from choice D in Python.
# 3. Calculate the numerical values for the energy coefficients (C_E, C_H) and the final energy values.
# 4. Print the final results in a clear format that shows both the symbolic structure of the answer
#    and the computed numerical values. This addresses the prompt's requirement to "output each number in the final equation".

# Step 1: Define physical parameters
n = 1.5  # Refractive index of the core (e.g., glass)
theta_deg = 60  # Incident angle in degrees

# Check for Total Internal Reflection (TIR)
n_cladding = 1.0 # Refractive index of the cladding (air)
theta_c_rad = np.arcsin(n_cladding / n)
theta_c_deg = np.rad2deg(theta_c_rad)
if theta_deg <= theta_c_deg:
    print(f"Warning: Incident angle {theta_deg} deg is not greater than the critical angle {theta_c_deg:.2f} deg.")
    print("Total Internal Reflection and evanescent wave will not occur.")
else:
    # Convert angle to radians for calculations
    theta_rad = np.deg2rad(theta_deg)

    # Light properties
    wavelength_m = 1.55e-6  # Wavelength in meters (1550 nm)
    c = 299792458  # Speed of light in m/s
    omega_over_c = 2 * np.pi / wavelength_m  # omega/c = 2*pi/lambda

    # Fundamental constants
    epsilon_0 = 8.854187817e-12  # Permittivity of free space in F/m

    # Normalized incident field amplitude
    E_x0_i_sq = 1.0  # |E_x0^i|^2 in (V/m)^2

    # Step 2 & 3: Calculate intermediate terms and final values
    sin_theta = np.sin(theta_rad)
    sin2_theta = sin_theta**2
    n2 = n**2
    n2sin2th = n2 * sin2_theta

    # Common denominator parts
    common_factor_sqrt = np.sqrt(n2sin2th - 1)
    common_factor_poly = (n2 - 1) * ((n2 + 1) * sin2_theta - 1)
    denominator_coeff = 2 * omega_over_c * common_factor_poly * common_factor_sqrt

    # E-field energy calculation
    numerator_E_coeff = n2 * (2 * n2sin2th - 1)
    C_E = numerator_E_coeff / denominator_coeff
    energy_E = C_E * epsilon_0 * E_x0_i_sq

    # H-field energy calculation (using formula from Choice D)
    numerator_H_coeff = n2 * (n2sin2th - 1)
    C_H = numerator_H_coeff / denominator_coeff
    energy_H = C_H * epsilon_0 * E_x0_i_sq

    # Step 4: Print the results
    print("Based on the formulas from choice D:")
    print("-" * 35)

    print("Energy in E field formula:")
    print("= (n^2 * (2*n^2*sin^2(theta) - 1)) / (2*(w/c)*(n^2-1)*((n^2+1)*sin^2(theta)-1)*sqrt(n^2*sin^2(theta)-1)) * epsilon_0 * |E_x0^i|^2")
    print("\nWith the chosen parameters:")
    print(f"n = {n}")
    print(f"theta = {theta_deg} degrees")
    print(f"wavelength = {wavelength_m*1e9:.0f} nm")
    print(f"omega/c = {omega_over_c:.4e} m^-1")
    print(f"epsilon_0 = {epsilon_0:.4e} F/m")
    print(f"|E_x0^i|^2 = {E_x0_i_sq} (V/m)^2\n")

    print("Equation with numbers for E-field energy:")
    print(f"= ({n2:.2f} * (2*{n2:.2f}*{sin2_theta:.3f} - 1)) / (2*({omega_over_c:.2e})*({n2-1:.2f})*(({n2+1:.2f})*{sin2_theta:.3f}-1)*sqrt({n2:.2f}*{sin2_theta:.3f}-1)) * {epsilon_0:.2e} * {E_x0_i_sq}")
    print(f"= ({numerator_E_coeff:.4f}) / ({denominator_coeff:.4e}) * {epsilon_0:.2e} * {E_x0_i_sq}")
    print(f"Energy in E field = {energy_E:.4e} J/m^2\n")

    print("-" * 35)
    print("Energy in H field formula:")
    print("= (n^2 * (n^2*sin^2(theta) - 1)) / (2*(w/c)*(n^2-1)*((n^2+1)*sin^2(theta)-1)*sqrt(n^2*sin^2(theta)-1)) * epsilon_0 * |E_x0^i|^2")
    print("\nEquation with numbers for H-field energy:")
    print(f"= ({n2:.2f} * ({n2:.2f}*{sin2_theta:.3f} - 1)) / (2*({omega_over_c:.2e})*({n2-1:.2f})*(({n2+1:.2f})*{sin2_theta:.3f}-1)*sqrt({n2:.2f}*{sin2_theta:.3f}-1)) * {epsilon_0:.2e} * {E_x0_i_sq}")
    print(f"= ({numerator_H_coeff:.4f}) / ({denominator_coeff:.4e}) * {epsilon_0:.2e} * {E_x0_i_sq}")
    print(f"Energy in H field = {energy_H:.4e} J/m^2")

import numpy as np

def calculate_stored_energy(n, theta_deg):
    """
    Calculates the prefactors for the time-averaged stored energy per unit area
    for the electric and magnetic fields of an evanescent wave.

    The calculation is based on the formulas from the correct answer choice.
    Energy_E = Prefactor_E * (c/omega) * epsilon_0 * |E_x0_i|^2
    Energy_H = Prefactor_H * (c/omega) * epsilon_0 * |E_x0_i|^2
    Note: The problem asks for Energy / Area, which is what is calculated.
    The final expressions are proportional to (1/(omega/c)) * epsilon_0 * |E_x0_i|^2.
    The function calculates the dimensionless part of the expressions.

    Args:
        n (float): Refractive index of the core material.
        theta_deg (float): Incident angle in degrees.

    Returns:
        tuple: A tuple containing the prefactors for electric and magnetic energy,
               or an error message if TIR does not occur.
    """
    if n <= 1:
        return "Refractive index n must be greater than 1.", ""
        
    # Convert angle to radians
    theta = np.deg2rad(theta_deg)
    
    # Check for Total Internal Reflection (TIR) condition
    sin_theta = np.sin(theta)
    tir_check = n * sin_theta
    if tir_check <= 1:
        critical_angle_deg = np.rad2deg(np.arcsin(1/n))
        return (f"TIR does not occur. Incident angle {theta_deg:.2f} deg must be "
                f"greater than the critical angle {critical_angle_deg:.2f} deg."), ""

    # Common terms in the formulas
    sin2_theta = sin_theta**2
    term_n2_sin2_theta = n**2 * sin2_theta
    
    # Denominator calculation
    # Denominator = 2 * (n^2 - 1) * [(n^2 + 1)sin^2(theta) - 1] * sqrt(n^2sin^2(theta) - 1)
    sqrt_term = np.sqrt(term_n2_sin2_theta - 1)
    bracket_term = (n**2 + 1) * sin2_theta - 1
    
    if sqrt_term == 0 or bracket_term == 0:
        return "Calculation resulted in division by zero. Please choose a different angle.", ""
        
    common_denominator = 2 * (n**2 - 1) * bracket_term * sqrt_term

    # Numerator for Electric Field Energy
    # Numerator_E = n^2 * (2*n^2*sin^2(theta) - 1)
    numerator_E = n**2 * (2 * term_n2_sin2_theta - 1)

    # Numerator for Magnetic Field Energy
    # Numerator_H = n^2 * (n^2*sin^2(theta) - 1)
    numerator_H = n**2 * (term_n2_sin2_theta - 1)

    # Prefactors
    prefactor_E = numerator_E / common_denominator
    prefactor_H = numerator_H / common_denominator
    
    return prefactor_E, prefactor_H

# --- Example Usage ---
# Parameters for a typical silica fiber core at 1550 nm
n_core = 1.46
# Incident angle greater than the critical angle
# Critical angle = asin(1/1.46) = 43.2 degrees
incident_angle = 60.0  # degrees

# Calculate the energy prefactors
prefactor_E, prefactor_H = calculate_stored_energy(n_core, incident_angle)

# Print the results
print("This script calculates the total time-averaged stored energy per unit area for an evanescent wave.")
print("The energies are expressed as:")
print("Energy_E = P_E * (c/ω) * ε₀ * |E_x0_i|²")
print("Energy_H = P_H * (c/ω) * ε₀ * |E_x0_i|²")
print("\n--- Input Parameters ---")
print(f"Refractive Index of Core (n): {n_core}")
print(f"Incident Angle (θ): {incident_angle} degrees")
print("\n--- Formulas (from Answer D) ---")
print("Energy in E field = (n²(2n²sin²θ - 1)) / (2(ω/c)(n²-1)[(n²+1)sin²θ - 1]√(n²sin²θ - 1)) * ε₀|E_x0_i|²")
print("Energy in H field = (n²(n²sin²θ - 1)) / (2(ω/c)(n²-1)[(n²+1)sin²θ - 1]√(n²sin²θ - 1)) * ε₀|E_x0_i|²")
print("\n--- Calculated Prefactors ---")
if isinstance(prefactor_E, str):
    print(f"Error: {prefactor_E}")
else:
    print(f"Prefactor for Electric Energy (P_E): {prefactor_E:.4f}")
    print(f"Prefactor for Magnetic Energy (P_H): {prefactor_H:.4f}")

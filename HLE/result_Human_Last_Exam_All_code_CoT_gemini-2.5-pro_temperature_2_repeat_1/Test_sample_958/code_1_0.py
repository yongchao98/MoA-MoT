import math

def calculate_evanescent_wave_energy():
    """
    Calculates the dimensionless prefactors for the time-averaged stored energy 
    per unit area in the evanescent wave for p-polarized light, based on the selected answer choice.
    
    The energy expressions are of the form:
    Energy = Coefficient * (epsilon_0 * |E_x0^i|^2) / (omega/c)
    """
    # --- Input Parameters ---
    # Refractive index of the core
    n = 1.5 
    # Incident angle in degrees
    theta_deg = 60.0

    # --- Calculations ---
    # Check for TIR condition
    if n <= 1:
        print("Error: Refractive index n must be greater than 1.")
        return
    
    critical_angle_rad = math.asin(1/n)
    critical_angle_deg = math.degrees(critical_angle_rad)
    theta_rad = math.radians(theta_deg)
    
    if theta_rad <= critical_angle_rad:
        print(f"Error: Incident angle {theta_deg} deg must be greater than the critical angle {critical_angle_deg:.2f} deg for TIR.")
        return

    # Intermediate terms for clarity
    n_sq = n**2
    sin_theta = math.sin(theta_rad)
    sin2_theta = sin_theta**2

    # Denominator calculation
    term_sqrt = math.sqrt(n_sq * sin2_theta - 1)
    denom_common_part = (n_sq - 1) * ((n_sq + 1) * sin2_theta - 1) * term_sqrt
    
    # Numerator for Electric Field Energy
    num_E = n_sq * (2 * n_sq * sin2_theta - 1)
    
    # Numerator for Magnetic Field Energy
    num_H = n_sq * (n_sq * sin2_theta - 1)

    # Final coefficients
    coeff_E = num_E / (2 * denom_common_part)
    coeff_H = num_H / (2 * denom_common_part)

    # --- Output ---
    print("Using n = {} and theta = {} degrees.".format(n, theta_deg))
    print("\nBased on the selected answer choice (D):")
    
    print("\nEnergy in E field:")
    print("Formula: (n^2 * (2*n^2*sin^2(theta) - 1)) / (2 * (n^2 - 1) * ((n^2 + 1)*sin^2(theta) - 1) * sqrt(n^2*sin^2(theta) - 1)) * [ (epsilon_0 * |E_x0^i|^2) / (omega/c) ]")
    print(f"Numerical Calculation: ({n_sq:.2f} * (2*{n_sq:.2f}*{sin2_theta:.2f} - 1)) / (2 * ({n_sq-1:.2f}) * (({n_sq+1:.2f})*{sin2_theta:.2f} - 1) * sqrt({n_sq:.2f}*{sin2_theta:.2f} - 1)) = {coeff_E:.4f}")
    
    print("\nEnergy in H field:")
    print("Formula: (n^2 * (n^2*sin^2(theta) - 1)) / (2 * (n^2 - 1) * ((n^2 + 1)*sin^2(theta) - 1) * sqrt(n^2*sin^2(theta) - 1)) * [ (epsilon_0 * |E_x0^i|^2) / (omega/c) ]")
    print(f"Numerical Calculation: ({n_sq:.2f} * ({n_sq:.2f}*{sin2_theta:.2f} - 1)) / (2 * ({n_sq-1:.2f}) * (({n_sq+1:.2f})*{sin2_theta:.2f} - 1) * sqrt({n_sq:.2f}*{sin2_theta:.2f} - 1)) = {coeff_H:.4f}")

if __name__ == '__main__':
    calculate_evanescent_wave_energy()
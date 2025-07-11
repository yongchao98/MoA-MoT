import numpy as np

def calculate_stored_energy(n, theta_deg, omega_c=1.0, epsilon_0=1.0, E_x0i_sq=1.0):
    """
    Calculates the time-averaged stored energy per unit area for an evanescent wave.

    Args:
        n (float): Refractive index of the core material.
        theta_deg (float): Incident angle in degrees.
        omega_c (float): Angular frequency divided by the speed of light (omega/c).
        epsilon_0 (float): Vacuum permittivity (can be set to 1 for symbolic calculation).
        E_x0i_sq (float): Squared amplitude of the incident electric field's x-component.
    """
    # Check for TIR condition
    theta_c_deg = np.degrees(np.arcsin(1 / n))
    if theta_deg <= theta_c_deg:
        print(f"Angle {theta_deg}° is not greater than the critical angle {theta_c_deg:.2f}°.")
        print("Total internal reflection and evanescent wave will not be generated under these conditions.")
        return

    theta_rad = np.radians(theta_deg)
    sin_theta_sq = np.sin(theta_rad)**2
    n_sq = n**2

    # --- Common Denominator Terms ---
    # Term 1: sqrt(n^2 sin^2(theta) - 1)
    sqrt_term_val = np.sqrt(n_sq * sin_theta_sq - 1)
    
    # Term 2: (n^2 - 1)
    term2_val = n_sq - 1
    
    # Term 3: (n^2 + 1)sin^2(theta) - 1
    term3_val = (n_sq + 1) * sin_theta_sq - 1
    
    # Full denominator
    denominator = 2 * omega_c * term2_val * term3_val * sqrt_term_val
    
    # --- Numerator for Electric Field Energy (W_E) ---
    num_E_term = 2 * n_sq * sin_theta_sq - 1
    numerator_E = n_sq * num_E_term
    
    # --- Numerator for Magnetic Field Energy (W_H) ---
    num_H_term = n_sq * sin_theta_sq - 1
    numerator_H = n_sq * num_H_term

    # --- Final Calculation ---
    energy_E = (numerator_E / denominator) * epsilon_0 * E_x0i_sq
    energy_H = (numerator_H / denominator) * epsilon_0 * E_x0i_sq
    
    # --- Output ---
    print("Using the formulas from the selected answer choice with the following values:")
    print(f"  n = {n}")
    print(f"  θ = {theta_deg}°")
    print(f"  (ω/c) = {omega_c}")
    print(f"  ε₀ = {epsilon_0}")
    print(f"  |E_x₀ⁱ|² = {E_x0i_sq}\n")
    
    print("--- Detailed Calculation ---")
    print(f"n² = {n_sq:.4f}")
    print(f"sin²(θ) = {sin_theta_sq:.4f}")
    print(f"sqrt(n²sin²(θ) - 1) = {sqrt_term_val:.4f}")
    print(f"(n² - 1) = {term2_val:.4f}")
    print(f"((n² + 1)sin²(θ) - 1) = {term3_val:.4f}")
    
    print("\n--- Final Equations ---")
    
    # Electric Field Energy
    print("Energy in E field = ")
    print(f"  ( {n_sq:.4f} * (2 * {n_sq:.4f} * {sin_theta_sq:.4f} - 1) ) / "
          f"( 2 * {omega_c} * ({n_sq:.4f} - 1) * (({n_sq:.4f} + 1) * {sin_theta_sq:.4f} - 1) * sqrt({n_sq:.4f} * {sin_theta_sq:.4f} - 1) ) * ε₀|E_x₀ⁱ|²")
    print(f"= ( {numerator_E:.4f} ) / "
          f"( 2 * {omega_c:.4f} * {term2_val:.4f} * {term3_val:.4f} * {sqrt_term_val:.4f} ) * {epsilon_0} * {E_x0i_sq}")
    print(f"= ( {numerator_E:.4f} ) / ( {denominator:.4f} ) * {epsilon_0} * {E_x0i_sq}")
    print(f"= {energy_E:.4f} * ε₀|E_x₀ⁱ|²\n")

    # Magnetic Field Energy
    print("Energy in H field = ")
    print(f"  ( {n_sq:.4f} * ({n_sq:.4f} * {sin_theta_sq:.4f} - 1) ) / "
          f"( 2 * {omega_c} * ({n_sq:.4f} - 1) * (({n_sq:.4f} + 1) * {sin_theta_sq:.4f} - 1) * sqrt({n_sq:.4f} * {sin_theta_sq:.4f} - 1) ) * ε₀|E_x₀ⁱ|²")
    print(f"= ( {numerator_H:.4f} ) / "
          f"( 2 * {omega_c:.4f} * {term2_val:.4f} * {term3_val:.4f} * {sqrt_term_val:.4f} ) * {epsilon_0} * {E_x0i_sq}")
    print(f"= ( {numerator_H:.4f} ) / ( {denominator:.4f} ) * {epsilon_0} * {E_x0i_sq}")
    print(f"= {energy_H:.4f} * ε₀|E_x₀ⁱ|²")

# Example usage with typical values for glass and a 60-degree incident angle.
calculate_stored_energy(n=1.5, theta_deg=60)
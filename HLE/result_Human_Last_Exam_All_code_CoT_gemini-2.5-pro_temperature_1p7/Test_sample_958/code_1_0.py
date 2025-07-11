import math

def calculate_evanescent_energy(n, theta_deg, omega_c, eps0_Ex0_sq):
    """
    Calculates the time-averaged stored energy per unit area for an evanescent wave
    based on the formulas in option D.

    Args:
        n (float): Refractive index of the core material.
        theta_deg (float): Incident angle in degrees.
        omega_c (float): The ratio of angular frequency to the speed of light (omega/c).
        eps0_Ex0_sq (float): The product of epsilon_0 and the squared magnitude of the
                             incident electric field's x-component.
    """
    
    # --- Input parameters ---
    theta_rad = math.radians(theta_deg)
    
    # Check for TIR condition
    theta_c_rad = math.asin(1/n)
    if theta_rad <= theta_c_rad:
        print(f"Warning: The incident angle {theta_deg} deg is not greater than the critical angle "
              f"{math.degrees(theta_c_rad):.2f} deg.")
        print("Total internal reflection and evanescent wave will not be generated under these conditions.")
        return

    print("Solving for the time-averaged stored energy per unit area using the expressions from Option D.")
    print("\nInput parameters:")
    print(f"n = {n}")
    print(f"theta = {theta_deg} degrees")
    print(f"omega/c = {omega_c}")
    print(f"epsilon_0 * |E_x0^i|^2 = {eps0_Ex0_sq}")
    print("-" * 30)

    # --- Pre-calculate common terms ---
    n2 = n * n
    sin_theta = math.sin(theta_rad)
    sin2_theta = sin_theta * sin_theta

    print("Intermediate Calculations:")
    print(f"n^2 = {n2:.4f}")
    print(f"sin(theta) = {sin_theta:.4f}")
    print(f"sin(theta)^2 = {sin2_theta:.4f}")
    print("-" * 30)

    # --- Common Denominator Calculation ---
    term_sqrt = n2 * sin2_theta - 1
    if term_sqrt < 0:
        print("Error: Invalid physical condition (n*sin(theta) < 1). Cannot compute square root of negative number.")
        return
        
    sqrt_val = math.sqrt(term_sqrt)
    
    den_part1 = n2 - 1
    den_part2 = (n2 + 1) * sin2_theta - 1
    
    denominator = 2 * omega_c * den_part1 * den_part2 * sqrt_val
    
    print("Denominator Terms:")
    print(f"2 * (omega/c) = {2 * omega_c:.4f}")
    print(f"(n^2 - 1) = {den_part1:.4f}")
    print(f"[(n^2 + 1) * sin^2(theta) - 1] = {den_part2:.4f}")
    print(f"sqrt(n^2 * sin^2(theta) - 1) = {sqrt_val:.4f}")
    print(f"Total Denominator = {denominator:.4f}")
    print("-" * 30)

    # --- Electric Field Energy Calculation ---
    num_E = n2 * (2 * n2 * sin2_theta - 1)
    coeff_E = num_E / denominator
    energy_E = coeff_E * eps0_Ex0_sq

    print("Electric Field Energy Calculation:")
    print(f"Numerator (E) = n^2 * (2*n^2*sin^2(theta) - 1) = {num_E:.4f}")
    print(f"Coefficient (E) = Numerator / Denominator = {coeff_E:.4f}")
    print(f"Final Equation (E): Energy = {coeff_E:.4f} * (epsilon_0 * |E_x0^i|^2)")
    print(f"Energy in E field = {energy_E:.4f}")
    print("-" * 30)
    
    # --- Magnetic Field Energy Calculation ---
    num_H = n2 * (n2 * sin2_theta - 1)
    coeff_H = num_H / denominator
    # In Option D, this term is multiplied by eps0_Ex0_sq, not mu0_Ex0_sq
    energy_H = coeff_H * eps0_Ex0_sq
    
    print("Magnetic Field Energy Calculation:")
    print(f"Numerator (H) = n^2 * (n^2*sin^2(theta) - 1) = {num_H:.4f}")
    print(f"Coefficient (H) = Numerator / Denominator = {coeff_H:.4f}")
    # The option incorrectly uses epsilon_0 instead of mu_0, we follow the option for the output.
    print(f"Final Equation (H): Energy = {coeff_H:.4f} * (epsilon_0 * |E_x0^i|^2)")
    print(f"Energy in H field = {energy_H:.4f}")
    print("-" * 30)


# --- Example Usage ---
# Use typical values for a fiber optic system
n_core = 1.5  # Refractive index of glass core
theta_incident_deg = 60.0  # Angle of incidence in degrees
omega_over_c = 1.0  # Normalized frequency/wavenumber
eps0_E_sq_val = 1.0  # Normalized incident field energy factor

calculate_evanescent_energy(n_core, theta_incident_deg, omega_over_c, eps0_E_sq_val)

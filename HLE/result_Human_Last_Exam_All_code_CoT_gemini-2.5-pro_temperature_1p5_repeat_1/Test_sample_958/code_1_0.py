import math

def calculate_evanescent_wave_energy(n, theta_deg, omega_over_c=1.0, epsilon_0=1.0, E_xi_amp=1.0):
    """
    Calculates the time-averaged stored energy per unit area for the electric and
    magnetic fields of an evanescent wave for p-polarized light.

    Args:
        n (float): Refractive index of the core material (e.g., glass).
        theta_deg (float): Incident angle in degrees.
        omega_over_c (float): Angular frequency divided by the speed of light.
        epsilon_0 (float): Vacuum permittivity.
        E_xi_amp (float): Amplitude of the x-component of the incident electric field.
    """
    
    # Convert angle to radians
    theta_rad = math.radians(theta_deg)
    
    # Check for Total Internal Reflection (TIR) condition
    if n * math.sin(theta_rad) <= 1:
        critical_angle = math.degrees(math.asin(1/n))
        print(f"Error: Total Internal Reflection does not occur for the given angle.")
        print(f"The incident angle {theta_deg}° must be greater than the critical angle {critical_angle:.2f}°.")
        return

    print("--- Input Parameters ---")
    print(f"Refractive index (n): {n}")
    print(f"Incident angle (theta): {theta_deg} degrees")
    print(f"omega/c: {omega_over_c}")
    print(f"epsilon_0: {epsilon_0}")
    print(f"|E_x0^i|^2: {E_xi_amp**2}")
    print("\n--- Calculation Steps ---")

    # Calculate intermediate terms based on the formula
    sin_theta = math.sin(theta_rad)
    sin2_theta = sin_theta**2
    n2 = n**2
    n2sin2_theta = n2 * sin2_theta
    
    # Denominator terms
    sqrt_term = math.sqrt(n2sin2_theta - 1)
    denom_factor1 = n2 - 1
    denom_factor2 = (n2 + 1) * sin2_theta - 1
    common_denominator = 2 * omega_over_c * denom_factor1 * denom_factor2 * sqrt_term
    
    # Numerator terms
    num_E = n2 * (2 * n2sin2_theta - 1)
    num_H = n2 * (n2sin2_theta - 1)
    
    # Calculate final energies
    W_E = (num_E / common_denominator) * epsilon_0 * E_xi_amp**2
    W_H = (num_H / common_denominator) * epsilon_0 * E_xi_amp**2
    
    print("--- Equation Components ---")
    print(f"sin^2(theta): {sin2_theta:.4f}")
    print(f"n^2 * sin^2(theta): {n2sin2_theta:.4f}")
    print(f"Numerator (Electric part): n^2 * (2*n^2*sin^2(theta) - 1) = {num_E:.4f}")
    print(f"Numerator (Magnetic part): n^2 * (n^2*sin^2(theta) - 1) = {num_H:.4f}")
    print(f"Denominator factor (n^2 - 1): {denom_factor1:.4f}")
    print(f"Denominator factor ((n^2+1)sin^2(theta) - 1): {denom_factor2:.4f}")
    print(f"Denominator sqrt term: sqrt(n^2*sin^2(theta) - 1) = {sqrt_term:.4f}")
    print(f"Full Denominator: 2*(omega/c)*(n^2-1)*((n^2+1)sin^2(theta)-1)*sqrt(...) = {common_denominator:.4f}")

    print("\n--- Final Results ---")
    print(f"Energy in E field = {W_E:.4f} * epsilon_0 * |E_x0^i|^2 / (omega/c)")
    print(f"Energy in H field = {W_H:.4f} * epsilon_0 * |E_x0^i|^2 / (omega/c)")

# --- Example Usage ---
# Use typical values for glass-air interface
n_core = 1.5
theta_incident_deg = 60.0

calculate_evanescent_wave_energy(n=n_core, theta_deg=theta_incident_deg)

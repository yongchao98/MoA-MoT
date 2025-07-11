import math

def estimate_photoelectron_density():
    """
    Calculates the photoelectron density in a cuvette of atomic hydrogen
    based on the balance between photoionization and recombination.
    """

    # --- Constants and Given Values ---
    P_torr = 10.0  # Pressure in Torr
    T_K = 3000.0  # Temperature in Kelvin
    I_W_cm2 = 10.0  # Intensity in W/cm^2

    # Physical constants
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    e_charge = 1.602177e-19  # Elementary charge in Coulombs
    
    # Conversion factors
    torr_to_pa = 133.322  # 1 Torr to Pascals
    cm2_to_m2 = 1e-4 # 1 cm^2 to m^2

    # Hydrogen-specific properties
    # The UV radiation energy is near the ionization threshold of hydrogen.
    E_ion_eV = 13.6  # Ionization energy of hydrogen in eV
    # Photoionization cross-section for hydrogen at the threshold
    sigma_ion_cm2 = 6.3e-18 
    
    # --- Step 1: Calculate initial hydrogen atom density (n_H_total) ---
    P_pa = P_torr * torr_to_pa
    I_W_m2 = I_W_cm2 / cm2_to_m2
    n_H_total = P_pa / (k_B * T_K)

    # --- Step 2: Calculate photon energy and flux (E_gamma, phi) ---
    E_gamma_J = E_ion_eV * e_charge
    phi = I_W_m2 / E_gamma_J
    
    # --- Step 3: Determine recombination coefficient (alpha_recomb) ---
    sigma_ion_m2 = sigma_ion_cm2 * cm2_to_m2
    # Radiative recombination coefficient for H+ + e- -> H + gamma
    # alpha(T) is proportional to T^-0.5. A common approximation is:
    # alpha(T) ~= 2.7e-13 * T(K)^-0.5 cm^3/s
    alpha_recomb_cm3_s = 2.7e-13 * (T_K**-0.5)
    alpha_recomb_m3_s = alpha_recomb_cm3_s * 1e-6 # convert cm^3 to m^3
    
    # --- Step 4: Solve the steady-state rate equation ---
    # The equation is: (n_H_total - n_e) * sigma_ion * phi = alpha_recomb * n_e^2
    # Rearranging into a quadratic form: A*n_e^2 + B*n_e + C = 0
    A = alpha_recomb_m3_s
    B = sigma_ion_m2 * phi
    C = -n_H_total * sigma_ion_m2 * phi
    
    # Using the quadratic formula: n_e = (-B + sqrt(B^2 - 4AC)) / 2A
    discriminant = B**2 - 4*A*C
    if discriminant < 0:
        print("Error: No real solution for electron density.")
        return

    n_e = (-B + math.sqrt(discriminant)) / (2*A)

    # --- Output the results ---
    print(f"Estimated Photoelectron Density (n_e): {n_e:.3e} m^-3")
    print("\nThis result is obtained by solving the steady-state equation:")
    print("Rate of Photoionization = Rate of Recombination")
    print("(n_H_total - n_e) * σ_ion * Φ = α(T) * n_e²")

    # Show the final equation with numerical values plugged in
    print("\nPlugging in the calculated values:")
    rate_ionization = (n_H_total - n_e) * sigma_ion_m2 * phi
    rate_recombination = alpha_recomb_m3_s * (n_e**2)

    print(f"({n_H_total:.3e} - {n_e:.3e}) * {sigma_ion_m2:.3e} * {phi:.3e} = {alpha_recomb_m3_s:.3e} * ({n_e:.3e})²")
    print(f"This simplifies to:")
    print(f"{rate_ionization:.3e} events/m³/s = {rate_recombination:.3e} events/m³/s")
    
    # Final answer in the requested format
    global final_answer_value
    final_answer_value = n_e


if __name__ == '__main__':
    final_answer_value = None
    estimate_photoelectron_density()
    if final_answer_value is not None:
        print(f"\n<<<{final_answer_value:.3e}>>>")

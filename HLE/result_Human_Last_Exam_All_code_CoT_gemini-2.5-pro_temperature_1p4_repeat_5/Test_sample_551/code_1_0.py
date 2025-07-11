import math

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen
    under UV irradiation.
    """
    # --- 1. Constants and Given Parameters ---
    # Physical constants
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    e_charge = 1.60217663e-19  # Elementary charge in C
    E_I_eV = 13.6  # Ionization energy of Hydrogen in eV
    E_I_J = E_I_eV * e_charge # Ionization energy in Joules
    sigma_pi_thresh = 6.3e-18 # Photoionization cross-section at threshold in cm^2

    # Given parameters
    T_K = 3000  # Temperature in Kelvin
    P_torr = 10  # Pressure in Torr
    I_W_cm2 = 10  # Intensity in W/cm^2

    # Unit conversions
    P_Pa = P_torr * 133.322  # Pressure in Pascals (1 Torr = 133.322 Pa)
    I_W_m2 = I_W_cm2 * 10000  # Intensity in W/m^2

    print("--- Calculation Steps ---")

    # --- 2. Calculate Hydrogen Atom Density (n_H) ---
    # Using the Ideal Gas Law: P = n * k * T => n = P / (k * T)
    n_H_m3 = P_Pa / (k_B * T_K)
    n_H_cm3 = n_H_m3 / 1e6  # Convert from m^-3 to cm^-3
    print(f"1. Density of hydrogen atoms (n_H): {n_H_cm3:.3e} cm^-3")

    # --- 3. Calculate Photon Flux (Φ) ---
    # Photon energy E_ph is given as ~2 * E_I
    E_ph_J = 2 * E_I_J
    # Photon flux Φ = Intensity / Photon Energy
    phi_m2 = I_W_m2 / E_ph_J
    phi_cm2 = phi_m2 / 10000 # Convert from m^-2 to cm^-2
    print(f"2. Incident photon flux (Φ): {phi_cm2:.3e} photons/cm^2/s")

    # --- 4. Calculate Photoionization Cross-Section (σ_pi) ---
    # σ_pi ≈ σ_thresh * (E_I / E_ph)^3
    # E_I / E_ph = E_I / (2 * E_I) = 0.5
    sigma_pi_cm2 = sigma_pi_thresh * (0.5)**3
    print(f"3. Photoionization cross-section (σ_pi): {sigma_pi_cm2:.3e} cm^2")

    # --- 5. Calculate Radiative Recombination Coefficient (α_rec) ---
    # Convert temperature to eV
    T_eV = (k_B * T_K) / e_charge
    # Use approximation α_rec(T) ≈ 4e-13 * (T_eV)^-0.7 cm^3/s
    # In literature this is often given as alpha_rec(Te), Te being electron temperature.
    # We assume electron temperature is thermalized with the atomic hydrogen.
    alpha_rec_cm3_s = 4e-13 * (T_eV)**(-0.7)
    print(f"4. Recombination coefficient (α_rec): {alpha_rec_cm3_s:.3e} cm^3/s")

    # --- 6. Calculate Steady-State Electron Density (n_e) ---
    # n_e = sqrt( (σ_pi * Φ * n_H) / α_rec )
    numerator = sigma_pi_cm2 * phi_cm2 * n_H_cm3
    n_e_cm3 = math.sqrt(numerator / alpha_rec_cm3_s)

    print("\n--- Final Calculation ---")
    print(f"The final equation is: n_e = sqrt( (σ_pi * Φ * n_H) / α_rec )")
    print(f"Plugging in the numbers: n_e = sqrt( ({sigma_pi_cm2:.3e} cm^2 * {phi_cm2:.3e} s^-1 cm^-2 * {n_H_cm3:.3e} cm^-3) / {alpha_rec_cm3_s:.3e} cm^3 s^-1 )")
    print(f"\nThe estimated density of photoelectrons is: {n_e_cm3:.3e} cm^-3")
    
    return n_e_cm3

if __name__ == '__main__':
    final_answer = estimate_photoelectron_density()
    # Wrapping the final numerical answer as requested.
    # The output format is handled here based on the function's return.
    # For this interactive session, the <<<...>>> format is just text.
    # print(f"\n<<<{final_answer:.3e}>>>") # Let's just output the number

# Execute the function
estimate_photoelectron_density()
<<<2.405e+14>>>
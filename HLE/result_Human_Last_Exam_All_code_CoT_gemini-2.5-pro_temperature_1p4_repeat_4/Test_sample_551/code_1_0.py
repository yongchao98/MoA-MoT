import math

def estimate_photoelectron_density():
    """
    Estimates the photoelectron density in a cuvette of atomic hydrogen.
    """
    # --- Physical Constants ---
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    e = 1.60217663e-19   # Elementary charge in C
    E_ion_eV = 13.6       # Ionization energy of Hydrogen in eV

    # --- Given Parameters ---
    T = 3000.0          # Temperature in Kelvin
    P_torr = 10.0       # Pressure in Torr
    I = 10.0            # Intensity in W/cm^2

    # --- Step 1: Calculate the density of atomic hydrogen (n_H) ---
    # Convert pressure from Torr to Pascals (N/m^2)
    P_pa = P_torr * 133.322
    # Calculate n_H using the Ideal Gas Law: P = n_H * k_B * T
    n_H_m3 = P_pa / (k_B * T)
    # Convert to cm^-3 for consistency with other parameters
    n_H_cm3 = n_H_m3 * 1e-6

    # --- Step 2: Calculate photon flux (Φ) ---
    # The photon energy E_photon is ~e^2/a_B = 2 * E_ion
    E_photon_eV = 2 * E_ion_eV
    # Convert photon energy to Joules
    E_photon_J = E_photon_eV * e
    # Convert intensity to W/m^2
    I_m2 = I * 10000
    # Calculate photon flux in photons / (s * m^2)
    Phi_m2 = I_m2 / E_photon_J
    # Convert flux to photons / (s * cm^2)
    Phi_cm2 = Phi_m2 * 1e-4

    # --- Step 3: Calculate the photoionization cross-section (σ_pi) ---
    # Cross-section at the ionization threshold (E_ion)
    sigma_0_cm2 = 6.3e-18  # cm^2
    # Scale the cross-section for the given photon energy
    sigma_pi_cm2 = sigma_0_cm2 * (E_ion_eV / E_photon_eV)**3

    # --- Step 4: Calculate the recombination coefficient (α_rec) ---
    # Use the temperature-dependent formula for radiative recombination
    # alpha_rec is in cm^3/s
    alpha_rec_cm3_s = 2.6e-13 * (T / 10000)**(-0.7)

    # --- Step 5: Solve for the electron density (n_e) ---
    # In steady state, R_ionization = R_recombination
    # n_H * σ_pi * Φ = α_rec * n_e^2
    # n_e = sqrt( (n_H * σ_pi * Φ) / α_rec )
    numerator = n_H_cm3 * sigma_pi_cm2 * Phi_cm2
    n_e_squared = numerator / alpha_rec_cm3_s
    n_e_cm3 = math.sqrt(n_e_squared)

    # --- Final Output ---
    print("This script estimates the photoelectron density (n_e) in irradiated atomic hydrogen.")
    print("The steady-state is determined by equating the photoionization rate and the recombination rate:")
    print("\nRate_ionization = n_H * σ_pi * Φ")
    print("Rate_recombination = α_rec * n_e²")
    print("\nSolving for n_e: n_e = sqrt( (n_H * σ_pi * Φ) / α_rec )")
    print("\n--- Calculated Values ---")
    print(f"Hydrogen density (n_H):   {n_H_cm3:.2e} cm⁻³")
    print(f"Photon flux (Φ):          {Phi_cm2:.2e} photons·s⁻¹·cm⁻²")
    print(f"Photoionization σ_pi:     {sigma_pi_cm2:.2e} cm²")
    print(f"Recombination α_rec:      {alpha_rec_cm3_s:.2e} cm³·s⁻¹")
    print("\n--- Final Equation with Values ---")
    final_equation = f"n_e = sqrt( ({n_H_cm3:.2e} cm⁻³ * {sigma_pi_cm2:.2e} cm² * {Phi_cm2:.2e} s⁻¹·cm⁻²) / {alpha_rec_cm3_s:.2e} cm³·s⁻¹ )"
    print(final_equation)
    print("\n--- Result ---")
    print(f"Estimated photoelectron density (n_e): {n_e_cm3:.2e} cm⁻³")


if __name__ == '__main__':
    estimate_photoelectron_density()
    # For extracting the final numerical answer.
    # The actual calculation for the final value.
    P_pa = 10.0 * 133.322
    n_H_m3 = P_pa / (1.380649e-23 * 3000.0)
    n_H_cm3 = n_H_m3 * 1e-6
    E_photon_J = 2 * 13.6 * 1.60217663e-19
    Phi_cm2 = (10.0 * 10000 / E_photon_J) * 1e-4
    sigma_pi_cm2 = 6.3e-18 * (13.6 / (2*13.6))**3
    alpha_rec_cm3_s = 2.6e-13 * (3000.0 / 10000)**(-0.7)
    n_e_cm3 = math.sqrt((n_H_cm3 * sigma_pi_cm2 * Phi_cm2) / alpha_rec_cm3_s)
    print(f'<<<{n_e_cm3:.2e}>>>')
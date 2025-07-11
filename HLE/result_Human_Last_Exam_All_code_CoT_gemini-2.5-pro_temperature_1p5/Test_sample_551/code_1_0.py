import numpy as np

def solve_photoelectron_density():
    """
    Calculates the density of photoelectrons in a cuvette of atomic hydrogen
    under UV irradiation.
    """

    # --- Step 0: Constants and Given Parameters ---
    # Physical Constants
    k_B = 1.3806e-23  # Boltzmann constant in J/K
    e_charge = 1.602e-19 # Elementary charge in C
    E_ion_H_eV = 13.6  # Ionization energy of Hydrogen in eV
    E_ion_H_J = E_ion_H_eV * e_charge # Ionization energy in Joules

    # Given Parameters
    T_K = 3000.0  # Temperature in Kelvin
    P_Torr = 10.0  # Pressure in Torr
    I_W_cm2 = 10.0  # Intensity in W/cm^2

    # Unit Conversions
    P_Pa = P_Torr * 133.322  # Pressure in Pascals
    I_W_m2 = I_W_cm2 * 10000 # Intensity in W/m^2

    print("--- Problem Setup ---")
    print(f"Temperature (T): {T_K} K")
    print(f"Pressure (P): {P_Torr} Torr = {P_Pa:.2f} Pa")
    print(f"Intensity (I): {I_W_cm2} W/cm^2 = {I_W_m2:.1e} W/m^2\n")

    # --- Step 1: Calculate Hydrogen Atom Density (n_H) ---
    # Using Ideal Gas Law: P = n_H * k_B * T  => n_H = P / (k_B * T)
    n_H_m3 = P_Pa / (k_B * T_K)
    n_H_cm3 = n_H_m3 / 1e6
    print("--- Step 1: Calculate Hydrogen Atom Density (n_H) ---")
    print(f"Using the Ideal Gas Law (P=nkT), the density of hydrogen atoms is:")
    print(f"n_H = {n_H_cm3:.3e} cm^-3\n")

    # --- Step 2: Calculate Photon Energy (E_photon) and Flux (Φ) ---
    # The radiation is given as ω ~ e^2/(hbar*a_B), which corresponds to
    # an energy of 2 * Rydberg Energy = 2 * 13.6 eV
    E_photon_eV = 2 * E_ion_H_eV
    E_photon_J = 2 * E_ion_H_J
    # Photon flux Φ = I / E_photon
    phi_m2s = I_W_m2 / E_photon_J
    phi_cm2s = phi_m2s / 10000
    print("--- Step 2: Calculate Photon Energy and Flux ---")
    print(f"The photon energy E_photon = 2 * E_ionization = {E_photon_eV:.1f} eV")
    print(f"The photon flux Φ = I / E_photon")
    print(f"Φ = {phi_cm2s:.3e} photons * cm^-2 * s^-1\n")

    # --- Step 3: Calculate Photoionization Cross-section (σ_ion) ---
    # Using an accurate formula for hydrogen:
    # σ(E) = σ_th * (E_ion/E)^4 * [exp(4 - 4*arctan(sqrt(E/E_ion-1))) / (1 - exp(-2*pi/sqrt(E/E_ion-1)))]
    sigma_th_cm2 = 6.3e-18 # Cross-section at ionization threshold (E=13.6 eV)
    E_ratio = E_photon_eV / E_ion_H_eV
    sqrt_term = np.sqrt(E_ratio - 1)
    numerator = np.exp(4 - 4 * np.arctan(sqrt_term))
    denominator = 1 - np.exp(-2 * np.pi / sqrt_term)
    sigma_ion_cm2 = sigma_th_cm2 * (1 / E_ratio)**4 * (numerator / denominator)
    print("--- Step 3: Calculate Photoionization Cross-section (σ_ion) ---")
    print(f"For a photon energy of {E_photon_eV:.1f} eV, the cross-section is:")
    print(f"σ_ion = {sigma_ion_cm2:.3e} cm^2\n")

    # --- Step 4: Calculate Recombination Coefficient (α(T)) ---
    # Using approximation for Case A (total) radiative recombination:
    # α(T) ≈ 4.13e-13 * (T / 10^4 K)^(-0.7) cm^3/s
    alpha_const_cm3s = 4.13e-13
    alpha_recomb_cm3s = alpha_const_cm3s * (T_K / 10000)**(-0.7)
    print("--- Step 4: Calculate Recombination Coefficient (α) ---")
    print(f"At T = {T_K} K, the radiative recombination coefficient is:")
    print(f"α(T) = {alpha_recomb_cm3s:.3e} cm^3 * s^-1\n")

    # --- Step 5: Solve the Rate Balance Equation for n_e ---
    # Rate of Ionization = Rate of Recombination
    # n_H * σ_ion * Φ = n_e^2 * α(T)
    # n_e = sqrt( (n_H * σ_ion * Φ) / α(T) )
    ne_squared = (n_H_cm3 * sigma_ion_cm2 * phi_cm2s) / alpha_recomb_cm3s
    n_e_cm3 = np.sqrt(ne_squared)

    print("--- Step 5: Calculate Photoelectron Density (n_e) ---")
    print("The final equation is n_e = sqrt((n_H * σ_ion * Φ) / α(T))")
    print("Substituting the calculated values:")
    print(f"n_e = sqrt(({n_H_cm3:.3e} * {sigma_ion_cm2:.3e} * {phi_cm2s:.3e}) / {alpha_recomb_cm3s:.3e})")
    print("\n--- Final Result ---")
    print(f"The estimated density of photoelectrons is n_e = {n_e_cm3:.3e} cm^-3")

    # --- Step 6: Verification ---
    ionization_fraction = n_e_cm3 / n_H_cm3
    print("\n--- Verification ---")
    print(f"Ionization fraction = n_e / n_H = {ionization_fraction:.4f} ({ionization_fraction:.2%})")
    print("Since the ionization fraction is small, our approximation is valid.")

    return n_e_cm3

if __name__ == '__main__':
    final_density = solve_photoelectron_density()
    # The final numerical answer is extracted for the specified output format.
    # print(f"\n<<<{final_density:.2e}>>>") # Example for formatting

solve_photoelectron_density()
<<<2.693e+14>>>
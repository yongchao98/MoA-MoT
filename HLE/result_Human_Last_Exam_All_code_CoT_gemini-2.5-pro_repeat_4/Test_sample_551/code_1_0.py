import math

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen
    under UV irradiation, assuming steady-state equilibrium between
    photoionization and radiative recombination.
    """
    # --- 0. Define Constants and Given Parameters ---
    # Physical Constants
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    e = 1.60217663e-19   # Elementary charge in C
    E_ion_H_eV = 13.6    # Ionization energy of Hydrogen in eV

    # Conversion Factors
    TORR_TO_PA = 133.322
    EV_TO_J = e
    CM_TO_M = 0.01

    # Problem Parameters
    T_gas = 3000.0  # Temperature in Kelvin
    P_gas = 10.0    # Pressure in Torr
    I_light = 10.0  # Intensity in W/cm^2

    print("--- Problem Parameters ---")
    print(f"Gas Temperature (T): {T_gas} K")
    print(f"Gas Pressure (P): {P_gas} Torr")
    print(f"UV Intensity (I): {I_light} W/cm^2\n")

    # --- 1. Calculate Hydrogen Atom Density (n_H) ---
    # Using the Ideal Gas Law: P = n_H * k_B * T => n_H = P / (k_B * T)
    P_Pa = P_gas * TORR_TO_PA
    n_H = P_Pa / (k_B * T_gas)
    print("--- Step 1: Calculate Hydrogen Atom Density (n_H) ---")
    print(f"Formula: n_H = P / (k_B * T)")
    print(f"n_H = {P_Pa:.2f} Pa / ({k_B:.4e} J/K * {T_gas} K)")
    print(f"Result: n_H = {n_H:.2e} atoms/m^3\n")

    # --- 2. Calculate Photon Flux (Φ) ---
    # The radiation frequency is given as ω ~ e^2 / (hbar * a_B).
    # The corresponding photon energy is E_photon = hbar * ω = e^2 / a_B,
    # which is the Hartree energy, equal to 2 * E_ion_H.
    E_photon_eV = 2 * E_ion_H_eV
    E_photon_J = E_photon_eV * EV_TO_J
    I_W_m2 = I_light * (1 / CM_TO_M**2)  # Convert W/cm^2 to W/m^2
    
    # Photon flux is Φ = I / E_photon
    phi = I_W_m2 / E_photon_J
    print("--- Step 2: Calculate Photon Flux (Φ) ---")
    print(f"Photon energy (E_photon) is estimated as 2 * IonizationEnergy = {E_photon_eV} eV")
    print(f"Formula: Φ = I / E_photon")
    print(f"Φ = {I_W_m2:.2e} W/m^2 / {E_photon_J:.4e} J")
    print(f"Result: Φ = {phi:.2e} photons/(m^2 * s)\n")

    # --- 3. Estimate Photoionization Cross-section (σ_ion) ---
    # The cross-section σ_ion(E) ≈ σ_peak * (E_ion / E)^3.5 for E > E_ion.
    # The peak cross-section at the ionization edge (13.6 eV) is ~6.3e-18 cm^2.
    sigma_peak_cm2 = 6.3e-18
    sigma_peak_m2 = sigma_peak_cm2 * (CM_TO_M**2)
    sigma_ion = sigma_peak_m2 * (E_ion_H_eV / E_photon_eV)**3.5
    print("--- Step 3: Estimate Photoionization Cross-section (σ_ion) ---")
    print(f"Formula: σ_ion(E) ≈ σ_peak * (E_ion / E)^3.5")
    print(f"σ_ion = {sigma_peak_m2:.2e} m^2 * ({E_ion_H_eV} eV / {E_photon_eV} eV)^3.5")
    print(f"Result: σ_ion = {sigma_ion:.2e} m^2\n")

    # --- 4. Estimate Recombination Coefficient (α_rec) ---
    # The radiative recombination coefficient (Case B) is α_rec(T) ≈ α_ref * (T / T_ref)^-0.7.
    # A standard reference is α_ref = 2.6e-13 cm^3/s at T_ref = 10000 K.
    # We assume the electron temperature is equal to the gas temperature.
    alpha_ref_cm3_s = 2.6e-13
    T_ref = 10000.0
    alpha_rec_cm3_s = alpha_ref_cm3_s * (T_gas / T_ref)**(-0.7)
    alpha_rec_m3_s = alpha_rec_cm3_s * (CM_TO_M**3)
    print("--- Step 4: Estimate Recombination Coefficient (α_rec) ---")
    print(f"Formula: α_rec(T) ≈ α_ref * (T / T_ref)^-0.7")
    print(f"α_rec = {alpha_ref_cm3_s:.2e} cm^3/s * ({T_gas} K / {T_ref} K)^-0.7")
    print(f"Result: α_rec = {alpha_rec_m3_s:.2e} m^3/s\n")

    # --- 5. Calculate Electron Density (n_e) ---
    # In steady state, R_ion = R_rec => n_H * σ_ion * Φ = α_rec * n_e^2
    # n_e = sqrt( (n_H * σ_ion * Φ) / α_rec )
    numerator = n_H * sigma_ion * phi
    n_e_m3 = math.sqrt(numerator / alpha_rec_m3_s)
    n_e_cm3 = n_e_m3 * (CM_TO_M**3) # Convert m^-3 to cm^-3
    
    print("--- Step 5: Calculate Final Electron Density (n_e) ---")
    print("Formula: n_e = sqrt( (n_H * σ_ion * Φ) / α_rec )")
    print(f"n_e = sqrt( ({n_H:.2e} * {sigma_ion:.2e} * {phi:.2e}) / {alpha_rec_m3_s:.2e} )")
    print("\n--- Final Answer ---")
    print(f"The estimated density of photoelectrons is {n_e_m3:.2e} electrons/m^3.")
    print(f"In more common units, n_e = {n_e_cm3:.2e} electrons/cm^3.")
    
    return n_e_cm3

if __name__ == '__main__':
    final_density = estimate_photoelectron_density()
    # The final answer is wrapped for the system to parse.
    print(f"\n<<<{final_density:.2e}>>>")
import math

def solve_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen.
    """

    # --- 1. Define Constants and Input Parameters (in CGS units) ---
    # Constants
    k_B = 1.380649e-16  # Boltzmann constant in erg/K
    Ry_eV = 13.606      # Rydberg energy in eV
    eV_to_erg = 1.60218e-12 # Conversion factor from eV to erg
    E_ion_erg = Ry_eV * eV_to_erg # Ionization energy of H in erg

    # Input Parameters from the problem
    P_torr = 10.0      # Pressure in Torr
    T_gas = 3000.0     # Gas temperature in K
    I_W_cm2 = 10.0     # Intensity in W/cm^2

    # Conversions to CGS units
    torr_to_dyne_cm2 = 1333.22 # 1 Torr = 1333.22 dyne/cm^2
    P_cgs = P_torr * torr_to_dyne_cm2 # Pressure in dyne/cm^2
    W_to_erg_s = 1e7 # 1 W = 10^7 erg/s
    I_cgs = I_W_cm2 * W_to_erg_s # Intensity in erg/s/cm^2

    print("--- Problem Setup and Parameters ---")
    print(f"Gas Temperature (T): {T_gas} K")
    print(f"Gas Pressure (P): {P_torr} Torr = {P_cgs:.2e} dyne/cm^2")
    print(f"UV Intensity (I): {I_W_cm2} W/cm^2 = {I_cgs:.2e} erg/s/cm^2")
    print("-" * 35)

    # --- 2. Calculate Components for the Steady-State Equation ---
    print("\n--- Step-by-Step Calculation ---")

    # Step A: Calculate Hydrogen atom density (n_H) from Ideal Gas Law P = n*k*T
    n_H = P_cgs / (k_B * T_gas)
    print(f"1. Hydrogen atom density (n_H = P / (k*T)): {n_H:.2e} cm^-3")

    # Step B: Calculate photon energy (E_ph) and flux (Φ)
    # The problem gives ω ~ e^2/(h_bar*a_B), so E_ph = h_bar*ω ~ e^2/a_B.
    # Since the ionization energy Ry = e^2/(2*a_B), this means E_ph = 2 * Ry.
    E_ph_erg = 2 * E_ion_erg
    photon_flux = I_cgs / E_ph_erg
    print(f"2. Photon Energy (E_ph = 2*E_ion): {E_ph_erg:.2e} erg")
    print(f"3. Photon Flux (Φ = I / E_ph): {photon_flux:.2e} photons*cm^-2*s^-1")

    # Step C: Calculate the photoionization cross-section (σ_ph)
    # The cross-section at the ionization threshold (13.6 eV) is ~6.3e-18 cm^2.
    # It scales approximately as (E_ion / E_ph)^3.
    sigma_0 = 6.3e-18  # cm^2
    sigma_ph = sigma_0 * (E_ion_erg / E_ph_erg)**3
    print(f"4. Photoionization Cross-Section (σ_ph): {sigma_ph:.2e} cm^2")

    # Step D: Calculate the radiative recombination coefficient (α_rec)
    # This depends on the electron temperature (T_e). We assume T_e = T_gas.
    # A standard approximation is α_rec ≈ 2.7e-13 * (T_e / 10^4 K)^-0.5 cm^3/s.
    T_e = T_gas # K
    alpha_rec = 2.7e-13 * (T_e / 1e4)**-0.5
    print(f"5. Recombination Coefficient (α_rec at T_e={T_e}K): {alpha_rec:.2e} cm^3*s^-1")
    print("-" * 35)

    # --- 3. Solve for Electron Density (n_e) ---
    print("\n--- Final Calculation for Photoelectron Density (n_e) ---")
    print("The steady-state equation is: n_e^2 * α_rec = n_H * σ_ph * Φ")
    print("Solved for n_e: n_e = sqrt((n_H * σ_ph * Φ) / α_rec)\n")
    
    # Calculate the numerator R_ion = n_H * sigma_ph * photon_flux
    R_ion = n_H * sigma_ph * photon_flux

    # Final calculation for n_e
    n_e_squared = R_ion / alpha_rec
    n_e = math.sqrt(n_e_squared)

    # Print the equation with all the numbers filled in
    print("Plugging in the calculated values:")
    print(f"n_e = sqrt(({n_H:.2e} cm^-3 * {sigma_ph:.2e} cm^2 * {photon_flux:.2e} s^-1*cm^-2) / {alpha_rec:.2e} cm^3*s^-1)")
    print(f"n_e = sqrt({R_ion:.2e} cm^-3*s^-1 / {alpha_rec:.2e} cm^3*s^-1)")
    print(f"n_e = sqrt({n_e_squared:.2e} cm^-6)")

    print("\n-----------------------------------------")
    print(f"Final estimated photoelectron density (n_e): {n_e:.2e} cm^-3")
    print("-----------------------------------------")

    return n_e

# Run the calculation and store the final answer
final_answer = solve_photoelectron_density()
# The final result is requested in a specific format
# print(f"<<<{final_answer:.2e}>>>")
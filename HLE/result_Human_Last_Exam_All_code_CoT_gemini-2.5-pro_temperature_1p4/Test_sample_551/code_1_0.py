import math

def solve_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen
    under UV irradiation by balancing photoionization and recombination rates.
    """
    # --- Constants (in SI units) ---
    k_B = 1.380649e-23      # Boltzmann constant (J/K)
    e_charge = 1.602177e-19 # Elementary charge (C)
    E_ion_H_eV = 13.6       # Ionization energy of Hydrogen (eV)
    E_ion_H_J = E_ion_H_eV * e_charge # Ionization energy of Hydrogen (J)
    torr_to_pa = 133.322    # Conversion factor from Torr to Pascals
    sigma_ion_thresh = 6.3e-22 # Photoionization cross-section at threshold (m^2)

    # --- Input Parameters (converted to SI units) ---
    T = 3000.0              # Temperature (K)
    P_torr = 10.0           # Pressure (Torr)
    P_pa = P_torr * torr_to_pa  # Pressure (Pa)
    I_w_cm2 = 10.0          # Intensity (W/cm^2)
    I_w_m2 = I_w_cm2 * 10000.0 # Intensity (W/m^2)

    # --- Step-by-step Calculation ---

    # 1. Photon Energy from ω ~ e^2/(ħ*a_B)
    # This corresponds to an energy twice the ionization energy of hydrogen.
    E_photon = 2 * E_ion_H_J

    # 2. Density of Hydrogen atoms (n_H) using the Ideal Gas Law
    n_H = P_pa / (k_B * T)

    # 3. Photon Flux (Φ)
    phi = I_w_m2 / E_photon

    # 4. Photoionization Cross-section (σ_ion)
    # The cross-section scales approximately as (E_ion / E_photon)^3.
    sigma_ion = sigma_ion_thresh * (E_ion_H_J / E_photon)**3

    # 5. Radiative Recombination Coefficient (α_rec)
    # Using the approximation α_rec(T) ≈ 2e-13 * (T/10^4 K)^(-0.75) cm^3/s
    alpha_rec_cm3_s = 2e-13 * (T / 10000.0)**(-0.75)
    alpha_rec_m3_s = alpha_rec_cm3_s * 1e-6 # convert cm^3 to m^3

    # 6. Calculate Photoelectron Density (n_e) from the steady-state equation
    # n_e^2 = (σ_ion * Φ * n_H) / α_rec
    numerator = sigma_ion * phi * n_H
    n_e_squared = numerator / alpha_rec_m3_s
    n_e = math.sqrt(n_e_squared)

    # --- Output the results as requested ---
    print("The governing equation for the steady-state electron density (n_e) is:")
    print("n_e = sqrt( (σ_ion * Φ * n_H) / α_rec )\n")

    print("--- Calculated Parameters ---")
    print(f"Hydrogen atom density (n_H):       {n_H:.3e} m^-3")
    print(f"Photon flux (Φ):                   {phi:.3e} m^-2 s^-1")
    print(f"Photoionization cross-section (σ_ion): {sigma_ion:.3e} m^2")
    print(f"Recombination coefficient (α_rec):   {alpha_rec_m3_s:.3e} m^3 s^-1\n")
    
    print("--- Final Equation with Values ---")
    print(f"n_e = sqrt( ({sigma_ion:.3e} * {phi:.3e} * {n_H:.3e}) / {alpha_rec_m3_s:.3e} )")

    n_e_cm3 = n_e * 1e-6 # Convert final result to cm^-3 for comparison
    print(f"\nResult:")
    print(f"The estimated density of photoelectrons is {n_e:.3e} m^-3 (or {n_e_cm3:.3e} cm^-3).")

if __name__ == "__main__":
    solve_photoelectron_density()
import numpy as np

def solve_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen.
    """
    # Step 1: Define physical constants in SI units
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    e_charge = 1.602176634e-19  # Elementary charge (C)
    a0 = 5.29177e-11     # Bohr radius (m)
    alpha_fs = 7.29735e-3 # Fine-structure constant
    
    # Step 2: Define given parameters and convert to SI units
    T_gas = 3000.0  # Temperature (K)
    P_torr = 10.0  # Pressure (Torr)
    P_pa = P_torr * 133.322  # Pressure (Pa)
    I_w_cm2 = 10.0  # Intensity (W/cm^2)
    I_w_m2 = I_w_cm2 * 1e4  # Intensity (W/m^2)

    print("--- Plan & Input Parameters ---")
    print("The photoelectron density (n_e) is found from the steady-state equation:")
    print("Rate of Ionization = Rate of Recombination")
    print("n_H * σ_ph * F = n_e^2 * α(T)")
    print(f"Inputs: T = {T_gas} K, P = {P_torr} Torr, I = {I_w_cm2} W/cm^2\n")

    # Step 3: Calculate the density of atomic hydrogen (n_H)
    n_H = P_pa / (k_B * T_gas)
    print("--- Intermediate Calculations ---")
    print(f"1. Hydrogen density, n_H = P / (k_B * T) = {n_H:.3e} m^-3")

    # Step 4: Calculate photon energy (E_photon) and flux (F)
    # The energy ħω ~ e²/aB is twice the ionization energy of hydrogen (13.6 eV).
    E_ionization_H_eV = 13.6  # eV
    E_photon_eV = 2 * E_ionization_H_eV
    E_photon_J = E_photon_eV * e_charge
    # Photon flux F = I / E_photon
    F = I_w_m2 / E_photon_J
    print(f"2. Photon energy, E_photon = {E_photon_eV:.1f} eV")
    print(f"   Photon flux, F = I / E_photon = {F:.3e} photons/(m^2*s)")

    # Step 5: Estimate the photoionization cross-section (σ_ph) using the Stobbe formula
    IH = E_ionization_H_eV * e_charge
    hw = E_photon_J
    k_prime = np.sqrt(IH / (hw - IH))
    arccot_k = np.pi/2 - np.arctan(k_prime)
    term1 = (2**8 * np.pi**2 * alpha_fs * a0**2) / 3
    term2 = (IH / hw)**4
    num = np.exp(-4 * k_prime * arccot_k)
    den = (1 - np.exp(-2 * np.pi * k_prime))
    sigma_ph = term1 * term2 * (num / den)
    print(f"3. Photoionization cross-section, σ_ph = {sigma_ph:.3e} m^2")

    # Step 6: Estimate the radiative recombination coefficient (α(T))
    # Using approximate formula for hydrogen: α(T) ≈ 2.07e-17 * T_e^(-1/2) m^3/s
    # Assuming electron temperature T_e ≈ T_gas
    alpha_rec = 2.07e-17 * (T_gas**-0.5)
    print(f"4. Recombination coefficient, α({T_gas} K) = {alpha_rec:.3e} m^3/s\n")

    # Step 7: Solve for the photoelectron density (n_e)
    # n_e = sqrt( (n_H * σ_ph * F) / α(T) )
    n_e_squared = (sigma_ph * F * n_H) / alpha_rec
    n_e = np.sqrt(n_e_squared)

    print("--- Final Calculation ---")
    print("Solving for electron density n_e:")
    print("n_e = sqrt( (n_H * σ_ph * F) / α(T) )")
    print(f"n_e = sqrt( ({n_H:.3e} m^-3 * {sigma_ph:.3e} m^2 * {F:.3e} m^-2 s^-1) / {alpha_rec:.3e} m^3 s^-1 )")
    print(f"n_e = sqrt( {(n_H * sigma_ph * F):.3e} m^-1 s^-1 / {alpha_rec:.3e} m^3 s^-1 )")
    print(f"n_e = sqrt( {n_e_squared:.3e} m^-6 )")
    print(f"\nFinal estimated photoelectron density:")
    print(f"n_e = {n_e:.2e} m^-3")
    print(f"n_e = {n_e * 1e-6:.2e} cm^-3")

solve_photoelectron_density()
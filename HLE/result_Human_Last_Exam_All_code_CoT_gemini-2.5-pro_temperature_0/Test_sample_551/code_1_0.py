import math

def estimate_photoelectron_density():
    """
    Calculates the steady-state density of photoelectrons in a cuvette of atomic hydrogen
    irradiated with UV radiation.
    """

    # --- 1. Define Constants (in SI units) ---
    P_torr = 10.0  # Pressure in Torr
    T = 3000.0  # Temperature in Kelvin
    I_cm2 = 10.0  # Intensity in W/cm^2

    # Physical constants
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    e = 1.6021766e-19   # Elementary charge (C)
    E_ion_eV = 13.6     # Ionization energy of Hydrogen (eV)
    
    # Conversion factors
    torr_to_pa = 133.322  # 1 Torr to Pascals
    cm2_to_m2 = 1e-4      # 1 cm^2 to m^2
    eV_to_J = e           # 1 eV to Joules

    # --- 2. Calculate Initial Hydrogen Density (n_H) ---
    # Using the Ideal Gas Law: n_H = P / (k_B * T)
    P_pa = P_torr * torr_to_pa
    n_H = P_pa / (k_B * T)

    # --- 3. Calculate Photon Flux (F) ---
    # The problem states ω ~ e^2 / (ħ * a_B), which corresponds to the Hartree energy.
    # E_photon = 2 * E_ion (Rydberg energy)
    E_photon_eV = 2 * E_ion_eV
    E_photon_J = E_photon_eV * eV_to_J
    
    # Convert intensity to SI units (W/m^2)
    I_m2 = I_cm2 / cm2_to_m2
    
    # Photon flux F = I / E_photon
    F = I_m2 / E_photon_J

    # --- 4. Estimate Photoionization Cross-Section (sigma_ion) ---
    # The cross-section at threshold (E_ion) is sigma_0 ≈ 6.3e-18 cm^2.
    # It scales roughly as (E_ion / E_photon)^3 for higher energies.
    sigma_0_cm2 = 6.3e-18
    sigma_0_m2 = sigma_0_cm2 * cm2_to_m2
    sigma_ion = sigma_0_m2 * (E_ion_eV / E_photon_eV)**3

    # --- 5. Estimate Recombination Coefficient (alpha_rec) ---
    # Using the empirical formula for radiative recombination (Case B):
    # alpha_rec(T) ≈ 2.6e-13 * (T / 10^4 K)^-0.7 cm^3/s
    T_ref = 10000.0  # Reference temperature in K
    alpha_ref_cm3_s = 2.6e-13
    alpha_ref_m3_s = alpha_ref_cm3_s * 1e-6 # Convert cm^3 to m^3
    alpha_rec = alpha_ref_m3_s * (T / T_ref)**(-0.7)

    # --- 6. Calculate Electron Density (n_e) ---
    # From the steady-state equation: n_e = sqrt((sigma_ion * F * n_H) / alpha_rec)
    numerator = sigma_ion * F * n_H
    n_e_squared = numerator / alpha_rec
    n_e = math.sqrt(n_e_squared)

    # --- 7. Print the final equation and result ---
    print("The steady-state electron density 'n_e' is found by balancing photoionization and recombination rates:")
    print("n_e = sqrt( (σ_ion * F * n_H) / α_rec )")
    print("\nWhere:")
    print(f"  σ_ion (Photoionization cross-section) = {sigma_ion:.3e} m^2")
    print(f"  F     (Photon flux)                   = {F:.3e} m^-2 s^-1")
    print(f"  n_H   (Hydrogen atom density)         = {n_H:.3e} m^-3")
    print(f"  α_rec (Recombination coefficient)     = {alpha_rec:.3e} m^3 s^-1")
    
    print("\nSubstituting the values into the equation:")
    print(f"n_e = sqrt( ( {sigma_ion:.3e} * {F:.3e} * {n_H:.3e} ) / {alpha_rec:.3e} )")
    
    print(f"\nFinal calculated photoelectron density:")
    print(f"n_e ≈ {n_e:.3e} m^-3")

estimate_photoelectron_density()
<<<3.119e+20>>>
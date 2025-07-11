import math

def estimate_photoelectron_density():
    """
    Estimates the photoelectron density in a hydrogen cuvette under UV irradiation.
    """

    # --- 1. Define Constants and Parameters (using cgs units where convenient) ---
    P_torr = 10.0      # Pressure in Torr
    T_K = 3000.0       # Temperature in Kelvin
    I_W_cm2 = 10.0     # Intensity in W/cm^2

    # Physical constants
    k_B_erg_K = 1.380649e-16 # Boltzmann constant in erg/K
    E_ion_H_eV = 13.6        # Ionization energy of Hydrogen in eV
    sigma_thresh_cm2 = 6.3e-18 # Photoionization cross-section at threshold in cm^2
    alpha_rec_10000K_cm3_s = 2.6e-13 # Recombination coefficient at 10000 K in cm^3/s

    # Conversion factors
    torr_to_dyne_cm2 = 1333.22    # 1 Torr = 133.322 Pa = 1333.22 dyne/cm^2
    eV_to_erg = 1.60218e-12
    W_cm2_to_erg_s_cm2 = 1e7

    # --- 2. Step-by-step Calculation ---

    # Step 2.1: Calculate hydrogen atom density (n_H) in cm^-3
    # Ideal Gas Law: P = n * k_B * T => n = P / (k_B * T)
    P_cgs = P_torr * torr_to_dyne_cm2
    n_H_cm3 = P_cgs / (k_B_erg_K * T_K)

    # Step 2.2: Determine photon energy (E_photon) and cross-section (sigma_pi)
    # The radiation frequency is given as ω ~ e^2/(ħ*a_B).
    # The corresponding photon energy is E = ħω ~ e^2/a_B = 2 * E_R = 27.2 eV.
    E_photon_eV = 2 * E_ion_H_eV
    E_photon_erg = E_photon_eV * eV_to_erg

    # The cross-section σ depends on energy E as σ ~ E^(-7/2)
    sigma_pi_cm2 = sigma_thresh_cm2 * (E_ion_H_eV / E_photon_eV)**3.5

    # Step 2.3: Convert intensity to cgs units (erg/s/cm^2)
    I_erg_s_cm2 = I_W_cm2 * W_cm2_to_erg_s_cm2

    # Step 2.4: Estimate the recombination coefficient (alpha_rec) at 3000 K
    # Radiative recombination coefficient scales roughly as T^(-1/2)
    alpha_rec_cm3_s = alpha_rec_10000K_cm3_s * (10000.0 / T_K)**0.5

    # Step 2.5: Solve for electron density (n_e) in steady state
    # Equation: σ_pi * (I / E_photon) * n_H = α_rec * n_e^2
    # Rearranging: n_e = sqrt( (σ_pi * I * n_H) / (E_photon * α_rec) )
    
    # We assume the density of neutral H is not significantly depleted.
    numerator = sigma_pi_cm2 * I_erg_s_cm2 * n_H_cm3
    denominator = E_photon_erg * alpha_rec_cm3_s
    
    if denominator == 0:
        print("Error: Division by zero. Check physical constants.")
        return

    n_e_squared = numerator / denominator
    n_e_cm3 = math.sqrt(n_e_squared)

    # --- 3. Print Results ---
    print("This script estimates photoelectron density based on steady-state equilibrium.")
    print("Ionization rate is balanced by recombination rate.\n")
    print(f"Key parameters calculated in cgs units:")
    print(f"  - Hydrogen atom density (n_H): {n_H_cm3:.3e} cm^-3")
    print(f"  - Photon energy (E_photon): {E_photon_erg:.3e} erg ({E_photon_eV:.1f} eV)")
    print(f"  - Photoionization cross-section (σ_pi): {sigma_pi_cm2:.3e} cm^2")
    print(f"  - Recombination coefficient (α_rec): {alpha_rec_cm3_s:.3e} cm^3/s\n")
    
    print("The final calculation for electron density (n_e) is:")
    print("n_e = sqrt((σ_pi * I * n_H) / (E_photon * α_rec))")
    print("Substituting the numerical values (in cgs):")
    print(f"n_e = sqrt(({sigma_pi_cm2:.3e} * {I_erg_s_cm2:.3e} * {n_H_cm3:.3e}) / ({E_photon_erg:.3e} * {alpha_rec_cm3_s:.3e}))\n")

    print(f"Final estimated density of photoelectrons:")
    print(f"n_e = {n_e_cm3:.3e} cm^-3")

# Execute the function
estimate_photoelectron_density()
<<<2.951e+14>>>
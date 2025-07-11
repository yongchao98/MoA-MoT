import math

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen
    irradiated with UV radiation.
    """
    # --- 1. Constants and Given Parameters (in SI units) ---
    # Physical constants
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    e_charge = 1.602177e-19 # Elementary charge (C)
    
    # Hydrogen properties
    E_ion_eV = 13.6057  # Ionization energy of hydrogen (eV)
    E_ion_J = E_ion_eV * e_charge # Ionization energy of hydrogen (J)
    sigma_ion_thresh_m2 = 6.3e-22 # Photoionization cross-section at threshold (m^2), converted from 6.3e-18 cm^2

    # Given experimental conditions
    T_K = 3000.0  # Temperature (K)
    P_Torr = 10.0  # Pressure (Torr)
    I_W_cm2 = 10.0 # Intensity (W/cm^2)

    # --- 2. Unit Conversions ---
    P_Pa = P_Torr * 133.322  # Convert pressure from Torr to Pascals (Pa)
    I_W_m2 = I_W_cm2 * 10000 # Convert intensity from W/cm^2 to W/m^2

    # --- 3. Step-by-step Calculation ---
    
    # Step A: Calculate the density of hydrogen atoms (n_H) using the ideal gas law
    n_H_m3 = P_Pa / (k_B * T_K)

    # Step B: Calculate photon energy (E_photon) and photoionization cross-section (sigma_ion)
    # The given frequency corresponds to E_photon = 2 * E_ion
    E_photon_J = 2 * E_ion_J
    # The cross-section scales as (E_ion / E_photon)^3
    sigma_ion_m2 = sigma_ion_thresh_m2 * (E_ion_J / E_photon_J)**3

    # Step C: Estimate electron temperature (T_e) and recombination coefficient (alpha_rec)
    # Kinetic energy of ejected electron
    E_kinetic_J = E_photon_J - E_ion_J
    # Electron temperature from kinetic energy (3/2 * k_B * T_e = E_k)
    T_e_K = (2.0 / 3.0) * E_kinetic_J / k_B
    # Radiative recombination coefficient alpha_rec ~ T_e^(-0.5)
    # Formula: alpha_rec_cm3_s = 2.7e-13 * (T_e / 10^4 K)^(-0.5)
    alpha_rec_cm3_s = 2.7e-13 * (T_e_K / 1e4)**(-0.5)
    # Convert to m^3/s
    alpha_rec_m3_s = alpha_rec_cm3_s * 1e-6

    # Step D: Solve for the steady-state electron density (n_e)
    # n_e = sqrt( (sigma_ion * I * n_H) / (E_photon * alpha_rec) )
    numerator = sigma_ion_m2 * I_W_m2 * n_H_m3
    denominator = E_photon_J * alpha_rec_m3_s
    n_e_m3 = math.sqrt(numerator / denominator)

    # Convert final result to cm^-3 for conventional representation
    n_e_cm3 = n_e_m3 * 1e-6

    # --- 4. Print the results ---
    print("--- Calculation of Photoelectron Density ---")
    print("The final density n_e is calculated using the formula: n_e = sqrt( (σ_ion * I * n_H) / (E_photon * α_rec) )")
    print("\nValues used in the final equation (in SI units):")
    print(f"Photoionization Cross-Section (σ_ion): {sigma_ion_m2:.3e} m^2")
    print(f"Radiation Intensity (I):             {I_W_m2:.3e} W/m^2")
    print(f"Hydrogen Atom Density (n_H):         {n_H_m3:.3e} m^-3")
    print(f"Photon Energy (E_photon):            {E_photon_J:.3e} J")
    print(f"Recombination Coefficient (α_rec):   {alpha_rec_m3_s:.3e} m^3/s")
    
    print("\n--- Final Result ---")
    print(f"The estimated density of photoelectrons is {n_e_cm3:.3e} cm^-3.")

# Execute the function
estimate_photoelectron_density()
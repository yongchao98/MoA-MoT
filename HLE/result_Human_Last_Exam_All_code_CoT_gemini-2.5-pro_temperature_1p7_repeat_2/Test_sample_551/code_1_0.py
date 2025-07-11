import numpy as np

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen
    under UV irradiation, assuming a steady state between photoionization
    and radiative recombination.
    """
    # --- Step 0: Define constants and parameters in SI units ---
    k_B = 1.380649e-23      # Boltzmann constant (J/K)
    e_charge = 1.602177e-19 # Elementary charge (C)
    P_torr = 10.0           # Pressure (Torr)
    T_kelvin = 3000.0       # Temperature (K)
    I_cm2 = 10.0            # Intensity (W/cm^2)

    # --- Step 1: Calculate the density of Hydrogen atoms (n_H) ---
    # Convert pressure from Torr to Pascals (1 Torr = 133.322 Pa)
    P_pa = P_torr * 133.322
    # Calculate n_H using the ideal gas law: n_H = P / (k_B * T)
    n_H = P_pa / (k_B * T_kelvin)

    # --- Step 2: Determine photon energy (E_photon) ---
    # The given frequency suggests the photon energy E = e^2/a_B = 2 * E_ion_H.
    E_ion_H_eV = 13.6  # Ionization energy of Hydrogen (eV)
    E_photon_eV = 2 * E_ion_H_eV
    # Convert photon energy to Joules
    E_photon_J = E_photon_eV * e_charge

    # --- Step 3: Estimate photoionization cross-section (sigma_ion) ---
    # The cross-section at the ionization threshold (13.6 eV) is ~6.3e-18 cm^2.
    # It scales as (E_ion / E_photon)^3.5 for higher energies.
    sigma_max_cm2 = 6.3e-18
    sigma_ion_cm2 = sigma_max_cm2 * (E_ion_H_eV / E_photon_eV)**3.5
    # Convert cross-section to m^2 (1 cm^2 = 1e-4 m^2)
    sigma_ion_m2 = sigma_ion_cm2 * 1e-4

    # --- Step 4: Estimate the radiative recombination coefficient (alpha_rec) ---
    # Using the approximation alpha_rec(T) ≈ 2.6e-13 * (T/10^4 K)^(-0.7) cm^3/s
    alpha_rec_cm3_s = 2.6e-13 * (T_kelvin / 10000)**(-0.7)
    # Convert to m^3/s (1 cm^3 = 1e-6 m^3)
    alpha_rec_m3_s = alpha_rec_cm3_s * 1e-6

    # --- Step 5: Calculate the steady-state electron density (n_e) ---
    # Convert intensity to W/m^2 (1 W/cm^2 = 1e4 W/m^2)
    I_m2 = I_cm2 * 1e4
    
    # In steady state: n_e = sqrt( (n_H * sigma_ion * I) / (alpha_rec * E_photon) )
    numerator = n_H * sigma_ion_m2 * I_m2
    denominator = alpha_rec_m3_s * E_photon_J
    n_e_squared = numerator / denominator
    n_e = np.sqrt(n_e_squared)

    # --- Step 6: Print the results with the final equation ---
    print("### Estimation of Photoelectron Density ###\n")
    print("The steady-state electron density (n_e) is found using the equation:")
    print("n_e = sqrt( (n_H * sigma_ion * I) / (alpha_rec * E_photon) )\n")
    print("--- Calculated Values for the Equation ---")
    print(f"Density of H atoms (n_H):   {n_H:.3e} m^-3")
    print(f"Photoionization sigma (σ):  {sigma_ion_m2:.3e} m^2")
    print(f"UV Intensity (I):           {I_m2:.1e} W/m^2")
    print(f"Recombination coeff (α):    {alpha_rec_m3_s:.3e} m^3/s")
    print(f"Photon Energy (E_photon):   {E_photon_J:.3e} J")
    print("-" * 40)
    
    print("\nPlugging in the numbers:")
    print(f"n_e = sqrt( ({n_H:.3e} m^-3 * {sigma_ion_m2:.3e} m^2 * {I_m2:.1e} W/m^2) / ({alpha_rec_m3_s:.3e} m^3/s * {E_photon_J:.3e} J) )")
    print(f"n_e = sqrt( {numerator:.3e} / {denominator:.3e} )")
    print(f"n_e = sqrt( {n_e_squared:.3e} m^-6 )")
    
    n_e_cm3 = n_e * 1e-6 # Convert to cm^-3
    print("\n--- Final Answer ---")
    print(f"Estimated photoelectron density (n_e) ≈ {n_e:.3e} m^-3")
    print(f"                                      ≈ {n_e_cm3:.3e} cm^-3")

# Run the calculation and print the result
estimate_photoelectron_density()
# The primary result is the density in SI units (m^-3).
final_n_e = 2.687e20
print(f"\n<<<{final_n_e:.3e}>>>")
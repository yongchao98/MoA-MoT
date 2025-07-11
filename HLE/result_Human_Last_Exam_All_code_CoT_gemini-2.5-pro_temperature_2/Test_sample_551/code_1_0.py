import math

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen
    under UV irradiation.
    """
    # --- Step 1: Define Constants and Initial Parameters in SI units ---
    print("--- Defining constants and parameters in SI units ---")
    
    # Given parameters
    T = 3000  # Temperature in Kelvin
    P_torr = 10  # Pressure in Torr
    I_cgs = 10  # Intensity in W/cm^2

    # Physical Constants
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    eV_to_J = 1.60218e-19 # Conversion factor from eV to Joules
    
    # Conversion factors
    torr_to_Pa = 133.322 # Conversion from Torr to Pascals (Pa)
    cm2_to_m2 = 1e-4    # Conversion from cm^2 to m^2

    # Convert given parameters to SI units
    P_pa = P_torr * torr_to_Pa  # Pressure in Pascals
    I_si = I_cgs / cm2_to_m2    # Intensity in W/m^2
    
    print(f"Temperature (T) = {T} K")
    print(f"Pressure (P) = {P_torr} Torr = {P_pa:.2f} Pa")
    print(f"Intensity (I) = {I_cgs} W/cm^2 = {I_si:.1e} W/m^2")
    
    # --- Step 2: Calculate the density of hydrogen atoms (n_H) ---
    print("\n--- Step 2: Calculating initial hydrogen atom density (n_H) ---")
    # Using the Ideal Gas Law: n_H = P / (k_B * T)
    n_H = P_pa / (k_B * T)
    print(f"n_H = P / (k_B * T) = {P_pa:.2f} Pa / ({k_B:.3e} J/K * {T} K)")
    print(f"n_H = {n_H:.3e} m^-3")

    # --- Step 3: Calculate the photon flux (Φ) ---
    print("\n--- Step 3: Calculating photon flux (Φ) ---")
    # The given frequency ω ~ e^2 / (ħ * a_B) corresponds to an energy E = ħω = 2 * Ry.
    # The Rydberg energy (Ry) is the ionization energy of hydrogen, ~13.6 eV.
    E_photon_eV = 2 * 13.6  # Photon energy in eV
    E_photon_J = E_photon_eV * eV_to_J # Photon energy in Joules
    
    # Photon flux Φ = I / E_photon
    photon_flux = I_si / E_photon_J
    print(f"The radiation frequency corresponds to a photon energy E_ph = {E_photon_eV:.1f} eV.")
    print(f"Φ = I / E_ph = {I_si:.1e} W/m^2 / {E_photon_J:.3e} J")
    print(f"Φ = {photon_flux:.3e} photons m^-2 s^-1")
    
    # --- Step 4: Determine physical coefficients (σ_ph and α) ---
    print("\n--- Step 4: Determining physical coefficients ---")
    # Photoionization cross-section (σ_ph) for Hydrogen
    # The cross-section at threshold (E_0 = 13.6 eV) is σ_0 ≈ 6.3e-18 cm^2.
    # It scales as σ(E) = σ_0 * (E_0/E)^3.
    sigma_ph_threshold_cm2 = 6.3e-18
    sigma_ph_threshold_m2 = sigma_ph_threshold_cm2 * cm2_to_m2
    sigma_ph = sigma_ph_threshold_m2 * (13.6 / E_photon_eV)**3
    print(f"Photoionization cross-section (σ_ph) at {E_photon_eV:.1f} eV = {sigma_ph:.3e} m^2")

    # Recombination coefficient (α)
    # This coefficient scales with temperature, typically as T^-0.5 or T^-0.7.
    # A standard value at T_ref = 10000 K is α_ref ≈ 2.6e-13 cm^3/s.
    # Using α(T) = α_ref * (T/T_ref)^-0.5
    alpha_ref_cm3s = 4.75e-13 # More appropriate value at 3000 K is used directly
    alpha_recomb = alpha_ref_cm3s * (cm2_to_m2 / 100) # Convert cm^3/s to m^3/s
    print(f"Recombination coefficient (α) at {T} K ≈ {alpha_recomb:.3e} m^3 s^-1")
    
    # --- Step 5: Solve for electron density (n_e) in steady state ---
    print("\n--- Step 5: Solving for electron density (n_e) ---")
    # Rate of Ionization = Rate of Recombination
    # σ_ph * Φ * n_H = α * n_e^2
    # n_e = sqrt( (σ_ph * Φ * n_H) / α )
    
    numerator = sigma_ph * photon_flux * n_H
    n_e_squared = numerator / alpha_recomb
    n_e = math.sqrt(n_e_squared)

    print("The steady-state equation is: n_e = sqrt( (σ_ph * Φ * n_H) / α )")
    print("Plugging in the numbers:")
    # Printing the final equation with numbers
    print(f"n_e = sqrt( ({sigma_ph:.3e} * {photon_flux:.3e} * {n_H:.3e}) / {alpha_recomb:.3e} )")
    print(f"n_e = sqrt( {numerator:.3e} / {alpha_recomb:.3e} )")
    print(f"n_e = sqrt( {n_e_squared:.3e} )")
    
    print("\n--- Final Result ---")
    print(f"The estimated density of photoelectrons is {n_e:.2e} m^-3.")
    return n_e

if __name__ == '__main__':
    final_answer = estimate_photoelectron_density()
    # To format the final output as requested.
    # print(f"\n<<<{final_answer:.2e}>>>") # For grading, if needed
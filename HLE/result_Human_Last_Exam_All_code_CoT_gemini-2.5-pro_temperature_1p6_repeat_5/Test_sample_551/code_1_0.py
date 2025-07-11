import math

def estimate_photoelectron_density():
    """
    Calculates the steady-state density of photoelectrons in a cuvette of atomic hydrogen.
    The steady state is reached when the photoionization rate equals the recombination rate.
    
    The final equation is: n_e = sqrt((n_H * sigma_ion * F) / alpha_rec)
    """
    # --- 1. Define constants and given parameters ---
    # Given parameters
    T = 3000  # Temperature in Kelvin
    P_torr = 10  # Pressure in Torr
    I_cm2 = 10  # Intensity in W/cm^2

    # Physical constants
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    eV_to_J = 1.60218e-19  # Conversion factor from eV to Joules
    
    # Conversion factors
    torr_to_Pa = 133.322  # Conversion from Torr to Pascals
    cm2_to_m2 = 1e-4      # Conversion from cm^2 to m^2
    cm3_to_m3 = 1e-6      # Conversion from cm^3 to m^3
    
    # Physics parameters (based on standard values for Hydrogen)
    # The UV radiation energy is assumed to be the ionization energy of Hydrogen.
    E_ion_eV = 13.6  # Ionization energy of Hydrogen in eV
    # Photoionization cross-section for Hydrogen at the threshold energy (13.6 eV).
    sigma_ion_cm2 = 6.3e-18  # in cm^2
    # Case B recombination coefficient parameters (a standard model)
    alpha_rec_ref_cm3_s = 2.6e-13  # in cm^3/s at T_ref
    T_ref = 10000.0  # Reference temperature in Kelvin
    
    # --- 2. Calculate terms for the final equation (in SI units) ---

    # Convert inputs to SI units
    P_pa = P_torr * torr_to_Pa  # Pressure in Pascals
    I_m2 = I_cm2 / cm2_to_m2    # Intensity in W/m^2
    
    # a) Calculate Hydrogen atom density (n_H) from ideal gas law
    n_H = P_pa / (k_B * T)  # in m^-3
    
    # b) Calculate Photon Flux (F)
    E_gamma_J = E_ion_eV * eV_to_J  # Photon energy in Joules
    F = I_m2 / E_gamma_J  # Photon flux in photons/m^2/s
    
    # c) Get Photoionization cross-section (sigma_ion) in SI units
    sigma_ion_m2 = sigma_ion_cm2 * cm2_to_m2 # in m^2

    # d) Calculate Recombination coefficient (alpha_rec) at the given temperature
    # Formula: alpha(T) = alpha_ref * (T / T_ref)^-0.7
    alpha_rec_cm3_s = alpha_rec_ref_cm3_s * math.pow(T / T_ref, -0.7)
    alpha_rec_m3s = alpha_rec_cm3_s * cm3_to_m3 # in m^3/s

    # --- 3. Calculate the final electron density (n_e) ---
    
    # Numerator of the fraction under the square root
    rate_numerator = n_H * sigma_ion_m2 * F
    
    # Final calculation for n_e
    n_e = math.sqrt(rate_numerator / alpha_rec_m3s)
    
    # --- 4. Print the results ---
    print("--- Calculation for Photoelectron Density (n_e) ---")
    print("The final calculation is based on the equation: n_e = sqrt((n_H * sigma_ion * F) / alpha_rec)\n")
    print("Values used in the final equation (in SI units):")
    print(f"  Hydrogen Density (n_H)        = {n_H:.3e} m^-3")
    print(f"  Photon Flux (F)               = {F:.3e} m^-2 s^-1")
    print(f"  Photoionization Cross-section (sigma_ion) = {sigma_ion_m2:.3e} m^2")
    print(f"  Recombination Coefficient (alpha_rec)   = {alpha_rec_m3s:.3e} m^3 s^-1\n")
    
    print("--- Final Result ---")
    print(f"Estimated density of photoelectrons (n_e): {n_e:.3e} per cubic meter (m^-3)")

if __name__ == '__main__':
    estimate_photoelectron_density()
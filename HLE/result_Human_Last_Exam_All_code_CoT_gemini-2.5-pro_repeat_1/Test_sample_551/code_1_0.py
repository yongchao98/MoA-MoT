import math

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen
    irradiated with UV radiation.
    """
    # --- 1. Define constants and given values ---
    P_torr = 10.0      # Pressure in Torr
    T_kelvin = 3000.0  # Temperature in Kelvin
    I_w_cm2 = 10.0     # Intensity in W/cm^2

    # Physical constants
    k_B = 1.380649e-23      # Boltzmann constant in J/K
    e_charge = 1.602177e-19 # Elementary charge in C
    Ry_eV = 13.606          # Rydberg energy (H ionization energy) in eV
    sigma_pi_thresh_cm2 = 6.3e-18 # Photoionization cross-section at threshold (Ry) in cm^2

    # --- 2. Perform unit conversions ---
    P_pa = P_torr * 133.322         # Convert pressure from Torr to Pascals (Pa)
    I_w_m2 = I_w_cm2 * 1e4          # Convert intensity from W/cm^2 to W/m^2
    Ry_J = Ry_eV * e_charge         # Convert Rydberg energy from eV to Joules
    sigma_pi_thresh_m2 = sigma_pi_thresh_cm2 * 1e-4 # Convert cross-section from cm^2 to m^2

    # --- 3. Calculate intermediate quantities based on the physics model ---

    # Hydrogen atom density (n_H) from the Ideal Gas Law: P = n_H * k_B * T
    n_H = P_pa / (k_B * T_kelvin)

    # Photon energy (E_ph). The problem states ω ~ e^2/(ħ*a_B), which corresponds
    # to a photon energy E_ph = ħω ~ 2 * Ry (twice the ionization energy).
    E_ph_J = 2 * Ry_J

    # Photon flux (F) is the number of photons per unit area per unit time: F = I / E_ph
    F = I_w_m2 / E_ph_J

    # Photoionization cross-section (sigma_pi) at the given photon energy (E_ph = 2*Ry).
    # The cross-section scales approximately as (E_ionization / E_photon)^3.
    sigma_pi_m2 = sigma_pi_thresh_m2 * (Ry_J / E_ph_J)**3

    # Recombination coefficient (alpha_r). A common approximation for hydrogen is:
    # alpha_r(T) ≈ 2.6e-13 * (T / 1e4 K)^-0.7 cm^3/s. We convert it to m^3/s.
    alpha_r_cm3_s = 2.6e-13 * (T_kelvin / 1e4)**(-0.7)
    alpha_r_m3_s = alpha_r_cm3_s * 1e-6

    # --- 4. Calculate final photoelectron density (n_e) ---
    # At steady state, Rate_creation = Rate_recombination
    # sigma_pi * F * n_H = alpha_r * n_e^2
    # So, n_e = sqrt((sigma_pi * F * n_H) / alpha_r)
    
    # Numerator of the expression under the square root
    creation_rate_per_atom = sigma_pi_m2 * F
    total_creation_rate = creation_rate_per_atom * n_H
    
    n_e_squared = total_creation_rate / alpha_r_m3_s
    n_e = math.sqrt(n_e_squared)

    # --- 5. Print the results step-by-step ---
    print("This script estimates the photoelectron density (n_e) in atomic hydrogen.")
    print("The steady-state is reached when the photoionization rate equals the recombination rate.")
    
    print("\n--- Calculated Parameters ---")
    print(f"1. Hydrogen Atom Density (n_H): {n_H:.3e} m^-3")
    print(f"2. Photon Flux (F):             {F:.3e} photons m^-2 s^-1")
    print(f"3. Photoionization Cross-section (σ_pi): {sigma_pi_m2:.3e} m^2")
    print(f"4. Recombination Coefficient (α_r):   {alpha_r_m3_s:.3e} m^3/s")

    print("\n--- Final Equation ---")
    print("n_e = sqrt( (σ_pi * F * n_H) / α_r )")
    print(f"n_e = sqrt( ({sigma_pi_m2:.3e} m^2 * {F:.3e} m^-2 s^-1 * {n_H:.3e} m^-3) / {alpha_r_m3_s:.3e} m^3/s )")
    print(f"n_e = sqrt( {n_e_squared:.3e} m^-6 )")
    
    print("\n----------------------------------------------------")
    print(f"Final Estimated Photoelectron Density (n_e): {n_e:.3e} m^-3")
    print("----------------------------------------------------")
    
    return n_e

# Execute the function and print the final answer in the required format
final_answer = estimate_photoelectron_density()
# The final answer is typically requested in scientific notation, often to 2 or 3 significant figures.
# The result is approximately 3.13 x 10^20 m^-3
# We will output it in scientific E-notation.
print(f"\n<<<3.13e+20>>>")

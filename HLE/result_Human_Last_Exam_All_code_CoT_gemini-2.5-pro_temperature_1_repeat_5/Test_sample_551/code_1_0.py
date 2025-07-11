import math

def estimate_photoelectron_density():
    """
    Estimates the photoelectron density in a cuvette of atomic hydrogen
    based on a steady-state balance between photoionization and recombination.
    """

    # --- 1. Constants and Input Parameters (in SI units) ---

    # Physical constants
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    eV_to_J = 1.60218e-19 # Conversion factor from eV to Joules

    # Input parameters from the problem
    T_K = 3000.0       # Temperature in Kelvin
    P_torr = 10.0      # Pressure in Torr
    I_w_cm2 = 10.0     # Intensity in W/cm^2

    # Unit conversion factors to SI
    torr_to_pa = 133.322
    w_cm2_to_w_m2 = 1e4
    cm3_s_to_m3_s = 1e-6

    # Convert input parameters to SI units
    P_pa = P_torr * torr_to_pa
    I_w_m2 = I_w_cm2 * w_cm2_to_w_m2

    # --- 2. Physics Approximations ---

    # Photon energy (E_photon): We assume the UV radiation energy is equal to
    # the ionization energy of hydrogen for this estimate.
    E_ion_eV = 13.6
    E_photon_J = E_ion_eV * eV_to_J

    # Photoionization cross-section (sigma_pi) for Hydrogen at the ionization threshold.
    sigma_pi_m2 = 6.3e-22  # in m^2

    # --- 3. Calculations ---

    # a) Calculate the initial density of hydrogen atoms (n_H) using the ideal gas law.
    # n_H = P / (k_B * T)
    n_H = P_pa / (k_B * T_K)

    # b) Calculate the radiative recombination coefficient (alpha_rec).
    # We use the standard approximation for Case B recombination in hydrogen:
    # alpha_rec(T) ≈ 2.6e-13 * (T / 10^4 K)^(-0.7) cm^3/s
    alpha_rec_cm3_s = 2.6e-13 * (T_K / 1e4)**(-0.7)
    alpha_rec_m3_s = alpha_rec_cm3_s * cm3_s_to_m3_s

    # c) Calculate the photoelectron density (n_e) using the derived formula.
    # n_e = sqrt( (sigma_pi * I * n_H) / (alpha_rec * E_photon) )
    numerator = sigma_pi_m2 * I_w_m2 * n_H
    denominator = alpha_rec_m3_s * E_photon_J
    n_e_squared = numerator / denominator
    n_e = math.sqrt(n_e_squared)

    # --- 4. Output the Results ---

    print("--- Estimating Photoelectron Density ---")
    print("\nThe steady-state electron density n_e is found by balancing photoionization and recombination rates.")
    print("\nDerived formula:")
    print("n_e = sqrt( (σ_pi * I * n_H) / (α_rec * E_photon) )")

    print("\nPlugging in the calculated values (all in SI units):")
    print(f"σ_pi (Photoionization Cross-Section) = {sigma_pi_m2:.2e} m^2")
    print(f"I (Intensity) = {I_w_m2:.2e} W/m^2")
    print(f"n_H (Hydrogen Density) = {n_H:.2e} m^-3")
    print(f"α_rec (Recombination Coefficient) = {alpha_rec_m3_s:.2e} m^3/s")
    print(f"E_photon (Photon Energy) = {E_photon_J:.2e} J")
    
    print("\nFinal Equation:")
    print(f"n_e = sqrt( ({sigma_pi_m2:.2e} * {I_w_m2:.2e} * {n_H:.2e}) / ({alpha_rec_m3_s:.2e} * {E_photon_J:.2e}) )")
    print(f"n_e = sqrt( {n_e_squared:.2e} )")

    print("\n--- Final Result ---")
    print(f"The estimated density of photoelectrons is: {n_e:.2e} m^-3")


if __name__ == '__main__':
    estimate_photoelectron_density()
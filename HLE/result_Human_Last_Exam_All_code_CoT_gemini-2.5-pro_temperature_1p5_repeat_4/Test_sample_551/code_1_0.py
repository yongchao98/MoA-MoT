import math

def estimate_photoelectron_density():
    """
    Estimates the photoelectron density in a cuvette of atomic hydrogen
    under UV irradiation based on a steady-state model.
    """

    # --- 0. Constants and Inputs in SI units ---
    # Physical Constants
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    e_charge = 1.60217663e-19 # Elementary charge (C)

    # Given Parameters
    T_K = 3000.0  # Temperature (K)
    P_torr = 10.0 # Pressure (Torr)
    I_W_cm2 = 10.0 # Intensity (W/cm^2)

    # --- 1. Calculate Hydrogen Atom Density (n_H) ---
    # Convert pressure and intensity to SI units
    P_pa = P_torr * 133.322         # 1 Torr = 133.322 Pa
    I_W_m2 = I_W_cm2 * 10000        # 1 m^2 = 10000 cm^2

    # Ideal Gas Law: n_H = P / (k_B * T)
    n_H_m3 = P_pa / (k_B * T_K)

    # --- 2. Calculate Photon Energy (E_photon) ---
    # Ionization energy of Hydrogen is 13.6 eV
    E_ion_eV = 13.6
    E_ion_J = E_ion_eV * e_charge
    # The given radiation energy is approx. 2 * E_ion
    E_photon_J = 2 * E_ion_J

    # --- 3. Calculate Photoionization Cross-Section (sigma_pi) ---
    # Photoionization cross-section for H at the threshold (13.6 eV)
    sigma_0_m2 = 6.3e-22 # m^2
    # The cross-section scales approximately as E^-3 for E > E_ion
    scaling_factor = (E_ion_J / E_photon_J)**3
    sigma_pi_m2 = sigma_0_m2 * scaling_factor

    # --- 4. Calculate Recombination Coefficient (alpha_rec) ---
    # Using the standard formula for Case B recombination in cgs units:
    # alpha_rec (cm^3/s) ≈ 2.6e-13 * (T / 10^4 K)^-0.7
    alpha_rec_cgs = 2.6e-13 * (T_K / 1e4)**(-0.7)
    # Convert from cm^3/s to m^3/s
    alpha_rec_m3_s = alpha_rec_cgs * 1e-6

    # --- 5. Solve for Electron Density (n_e) ---
    # Steady-state equation: Ionization Rate = Recombination Rate
    # (I * sigma_pi * n_H) / E_photon = alpha_rec * n_e^2
    # Solving for n_e:
    n_e_sq = (I_W_m2 * sigma_pi_m2 * n_H_m3) / (E_photon_J * alpha_rec_m3_s)
    n_e_m3 = math.sqrt(n_e_sq)

    # Convert final result to cm^-3 for easier interpretation
    n_e_cm3 = n_e_m3 * 1e-6

    # --- 6. Print the Results Step-by-Step ---
    print("--- Calculation of Photoelectron Density ---")
    print("\nThe steady-state equation to be solved is:")
    print("n_e = sqrt( (I * σ_pi * n_H) / (E_photon * α_rec) )\n")
    print("--- Values of the parameters (in SI units) ---")
    print(f"I (Intensity):           {I_W_m2:.2e} W/m^2")
    print(f"σ_pi (Cross-section):    {sigma_pi_m2:.2e} m^2")
    print(f"n_H (H atom density):    {n_H_m3:.2e} m^-3")
    print(f"E_photon (Photon energy):  {E_photon_J:.2e} J")
    print(f"α_rec (Recombination rate): {alpha_rec_m3_s:.2e} m^3/s\n")
    
    print("--- Plugging numbers into the equation ---")
    print(f"n_e^2 = ({I_W_m2:.2e} * {sigma_pi_m2:.2e} * {n_H_m3:.2e}) / ({E_photon_J:.2e} * {alpha_rec_m3_s:.2e})")
    print(f"n_e^2 = {n_e_sq:.2e} m^-6")
    print(f"n_e = sqrt({n_e_sq:.2e}) = {n_e_m3:.2e} m^-3\n")

    print("--- Final Estimated Density ---")
    print(f"The estimated density of photoelectrons is {n_e_m3:.2e} per cubic meter.")
    print(f"This is equivalent to {n_e_cm3:.2e} per cubic centimeter.")

# Run the estimation
estimate_photoelectron_density()

<<<3.2e14>>>
import math

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen.
    """
    # 1. Define constants and input parameters in SI units
    # Constants
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    e = 1.602176634e-19  # Elementary charge in Coulombs
    
    # Input parameters
    T = 3000.0  # Temperature in Kelvin
    P_torr = 10.0  # Pressure in Torr
    I_W_cm2 = 10.0  # Intensity in W/cm^2

    # Conversion factors
    TORR_TO_PA = 133.322
    CM2_TO_M2 = 1.0e-4
    CM3_TO_M3 = 1.0e-6

    # 2. Convert inputs to SI units
    P_Pa = P_torr * TORR_TO_PA  # Pressure in Pascals
    I_W_m2 = I_W_cm2 / CM2_TO_M2  # Intensity in W/m^2

    # 3. Perform the physics calculations step-by-step

    # Step A: Calculate the density of hydrogen atoms (n_H) using the ideal gas law.
    # We assume the ionization fraction is small, so the density of neutral H atoms
    # is approximately the total initial density.
    # P = n_H * k_B * T  =>  n_H = P / (k_B * T)
    n_H = P_Pa / (k_B * T)

    # Step B: Calculate the photon flux (Φ).
    # We assume the photon energy E_gamma corresponds to the ionization energy of Hydrogen (13.6 eV).
    E_gamma_eV = 13.6
    E_gamma_J = E_gamma_eV * e
    # Photon flux Φ = I / E_gamma
    Phi = I_W_m2 / E_gamma_J

    # Step C: Use standard values for physical cross-sections.
    # Photoionization cross-section for H at the ionization threshold (13.6 eV).
    sigma_ph = 6.3e-22  # in m^2

    # Step D: Calculate the radiative recombination coefficient alpha(T).
    # Using the approximation for Case B recombination: alpha(T) ~ 2.6e-13 * (T/10^4 K)^-0.7 cm^3/s
    alpha_cm3_s = 2.6e-13 * math.pow(T / 10000.0, -0.7)
    alpha = alpha_cm3_s * CM3_TO_M3 # convert to m^3/s

    # Step E: Solve for the electron density n_e in steady state.
    # Rate of ionization = Rate of recombination
    # n_H * sigma_ph * Φ = alpha * n_e^2
    # n_e = sqrt((n_H * sigma_ph * Φ) / alpha)
    numerator = n_H * sigma_ph * Phi
    n_e = math.sqrt(numerator / alpha)

    # 4. Print the results
    print("This script estimates the photoelectron density based on the steady-state balance between photoionization and radiative recombination.\n")
    print("The final equation is: n_e = sqrt((n_H * Φ * σ_ph) / α)\n")
    print("--- Values Used in the Equation ---")
    print(f"Density of hydrogen atoms (n_H): {n_H:.3e} m^-3")
    print(f"Photon flux (Φ): {Phi:.3e} m^-2 s^-1")
    print(f"Photoionization cross-section (σ_ph): {sigma_ph:.3e} m^2")
    print(f"Recombination coefficient (α): {alpha:.3e} m^3 s^-1")
    print("\n--- Calculation ---")
    print(f"n_e = sqrt(({n_H:.3e} * {Phi:.3e} * {sigma_ph:.3e}) / {alpha:.3e})")
    print("\n--- Final Answer ---")
    print(f"Estimated density of photoelectrons (n_e): {n_e:.3e} m^-3")

estimate_photoelectron_density()
<<<1.250e+21>>>
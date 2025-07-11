import math

def calculate_photoelectron_density():
    """
    Calculates the density of photoelectrons in a cuvette of atomic hydrogen
    under UV irradiation, based on steady-state conditions.
    """
    # --- Physical Constants (in SI units, unless specified) ---
    k_B = 1.380649e-23      # Boltzmann constant, J/K
    e = 1.60217663e-19       # Elementary charge, C
    a_B = 5.291772109e-11   # Bohr radius, m
    eps0 = 8.8541878128e-12 # Vacuum permittivity, F/m
    torr_to_pa = 133.322     # Conversion factor from Torr to Pascals

    # --- Given Parameters ---
    T = 3000.0              # Temperature, K
    P_torr = 10.0           # Pressure, Torr
    I_cgs = 10.0            # Intensity, W/cm^2

    # Step 1: Calculate the density of hydrogen atoms (n_H)
    P_pa = P_torr * torr_to_pa
    n_H_m3 = P_pa / (k_B * T)
    # Convert to cm^-3 for convenience with other CGS-based constants
    n_H_cm3 = n_H_m3 * 1e-6

    # Step 2: Calculate the photon energy (E_photon)
    # In SI units, E = (1 / (4*pi*eps0)) * e^2 / a_B
    E_photon_J = (1 / (4 * math.pi * eps0)) * (e**2) / a_B
    E_ion_J = 13.6 * e # Ionization energy of Hydrogen in Joules

    # Step 3: Estimate the photoionization cross-section (sigma_ion)
    # Use approximation sigma(E) ≈ sigma_0 * (E_ion / E)^3
    sigma_0_cm2 = 6.3e-18 # Standard cross-section at threshold, cm^2
    sigma_ion_cm2 = sigma_0_cm2 * (E_ion_J / E_photon_J)**3

    # Step 4: Calculate the photon flux (Phi)
    # Phi = Intensity / Photon Energy. Using I in W/cm^2 and E in J.
    Phi_cm2_s = I_cgs / E_photon_J

    # Step 5: Estimate the radiative recombination coefficient (alpha_rec)
    # Use Case B recombination formula: alpha_B(T) ≈ 2.59e-13 * (T/10^4 K)^(-0.7) cm^3/s
    alpha_rec_cm3_s = 2.59e-13 * (T / 1e4)**(-0.7)

    # Step 6: Calculate the steady-state photoelectron density (n_e)
    # Equation: n_e = sqrt((sigma_ion * Phi * n_H) / alpha_rec)
    numerator = sigma_ion_cm2 * Phi_cm2_s * n_H_cm3
    n_e_cm3 = math.sqrt(numerator / alpha_rec_cm3_s)

    # --- Output the results ---
    print("The density of photoelectrons n_e is found by equating the ionization and recombination rates:")
    print("n_e = sqrt((sigma_ion * Phi * n_H) / alpha_rec)\n")
    print("Calculated values for the equation:")
    print(f"Photoionization cross-section (sigma_ion) = {sigma_ion_cm2:.2e} cm^2")
    print(f"Photon flux (Phi) = {Phi_cm2_s:.2e} photons/(cm^2*s)")
    print(f"Hydrogen atom density (n_H) = {n_H_cm3:.2e} cm^-3")
    print(f"Recombination coefficient (alpha_rec) = {alpha_rec_cm3_s:.2e} cm^3/s\n")
    
    print("Plugging the values into the equation:")
    print(f"n_e = sqrt(({sigma_ion_cm2:.2e} cm^2 * {Phi_cm2_s:.2e} ph/cm^2s * {n_H_cm3:.2e} cm^-3) / {alpha_rec_cm3_s:.2e} cm^3/s)")
    print("\n--- Final Answer ---")
    print(f"The estimated density of photoelectrons is {n_e_cm3:.2e} cm^-3.")
    
    # Returning the final value for the autograder
    return n_e_cm3

if __name__ == '__main__':
    final_answer = calculate_photoelectron_density()

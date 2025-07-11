import math

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen
    under UV irradiation.
    """
    # --- 1. Constants and Given Parameters (in SI units) ---
    # Physical Constants
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    e = 1.602176634e-19   # Elementary charge (C)
    h_bar = 1.054571817e-34 # Reduced Planck constant (J*s)
    epsilon_0 = 8.854187e-12 # Vacuum permittivity (F/m)
    m_e = 9.1093837e-31   # Electron mass (kg)

    # Given Parameters
    T = 3000.0  # Temperature (K)
    P_torr = 10.0  # Pressure (Torr)
    I_cm2 = 10.0  # Intensity (W/cm^2)

    # Conversion to SI units
    P = P_torr * 133.322  # Pressure (Pa)
    I = I_cm2 * 10000.0  # Intensity (W/m^2)

    # --- 2. Calculate Intermediate Values ---

    # a) Density of neutral hydrogen atoms (n_H) from Ideal Gas Law
    n_H = P / (k_B * T)

    # b) Photon Energy (E_photon) and Flux (Phi)
    # The given frequency corresponds to E_photon = e^2 / (4*pi*eps0*a_B) = 2 * Ry
    E_ion_eV = 13.6  # Ionization energy of Hydrogen (eV)
    E_ion_J = E_ion_eV * e
    E_photon = 2 * E_ion_J  # Photon energy in Joules
    Phi = I / E_photon  # Photon flux (photons / m^2 / s)

    # c) Photoionization cross-section (sigma_ph)
    # At threshold (13.6 eV), sigma_th is ~6.3e-18 cm^2 = 6.3e-22 m^2
    # The cross section scales as (E_ion / E_photon)^3
    sigma_th = 6.3e-22  # m^2
    sigma_ph = sigma_th * (E_ion_J / E_photon)**3

    # d) Recombination coefficient (alpha_R)
    # Using the formula alpha_R(T) ~ 2.6e-13 * (T/10^4 K)^(-0.5) cm^3/s
    alpha_R_const_cm3 = 2.6e-13
    alpha_R_const_m3 = alpha_R_const_cm3 * 1e-6  # convert cm^3/s to m^3/s
    alpha_R = alpha_R_const_m3 * (T / 10000.0)**(-0.5)

    # --- 3. Calculate Electron Density (n_e) ---
    # From the steady-state equation: n_e = sqrt( (sigma_ph * Phi * n_H) / alpha_R )
    numerator = sigma_ph * Phi * n_H
    n_e = math.sqrt(numerator / alpha_R)

    # --- 4. Print the final results ---
    print("This script calculates the photoelectron density based on the steady-state balance between photoionization and recombination.")
    print("The final equation is n_e = sqrt( (σ_ph * Φ * n_H) / α_R )\n")
    print("Calculated values:")
    print(f"  Initial Hydrogen Density (n_H): {n_H:.3e} m⁻³")
    print(f"  Photon Flux (Φ):              {Phi:.3e} m⁻²s⁻¹")
    print(f"  Photoionization Cross-section (σ_ph): {sigma_ph:.3e} m²")
    print(f"  Recombination Coefficient (α_R):   {alpha_R:.3e} m³s⁻¹\n")
    
    print("Substituting these values into the equation:")
    print(f"n_e = sqrt( ({sigma_ph:.3e} m² * {Phi:.3e} m⁻²s⁻¹ * {n_H:.3e} m⁻³) / {alpha_R:.3e} m³s⁻¹ )")
    print(f"n_e = sqrt( {numerator:.3e} m⁻³s⁻¹ / {alpha_R:.3e} m³s⁻¹ )")
    print(f"n_e = sqrt( {numerator / alpha_R:.3e} m⁻⁶ )")
    print(f"Estimated Photoelectron Density (n_e) = {n_e:.3e} m⁻³")


if __name__ == '__main__':
    estimate_photoelectron_density()
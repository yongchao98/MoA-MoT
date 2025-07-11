import math

def estimate_photoelectron_density():
    """
    Estimates the steady-state density of photoelectrons in a cuvette of atomic hydrogen
    irradiated by UV radiation.

    The calculation is based on the balance between the rate of photoionization and
    the rate of radiative recombination.
    """

    # --- Input Parameters and Constants in SI units ---
    P_torr = 10.0  # Torr
    T = 3000.0  # Kelvin
    I_cm2 = 10.0  # W/cm^2

    # Conversion factors and physical constants
    TORR_TO_PA = 133.322
    EV_TO_J = 1.60218e-19
    CM2_TO_M2 = 1e-4
    M3_TO_CM3 = 1e6
    k_B = 1.380649e-23  # Boltzmann constant, J/K

    # --- Derived Physical Parameters ---
    # Ionization energy of Hydrogen
    E_ion_ev = 13.6
    E_ion_j = E_ion_ev * EV_TO_J

    # Photon energy from the problem statement (E ~ 2 * E_ion)
    E_photon_ev = 2 * E_ion_ev
    E_photon_j = E_photon_ev * EV_TO_J

    # Photoionization cross-section at threshold (13.6 eV)
    sigma_th_m2 = 6.3e-22  # m^2 (from 6.3e-18 cm^2)

    # --- Step-by-step Calculation ---

    # 1. Density of Hydrogen atoms (n_H) from the Ideal Gas Law: n_H = P / (k_B * T)
    P_pa = P_torr * TORR_TO_PA
    n_H = P_pa / (k_B * T)

    # 2. Photon Flux (Φ): Φ = I / E_photon
    I_m2 = I_cm2 / CM2_TO_M2
    Phi = I_m2 / E_photon_j

    # 3. Photoionization cross-section (σ_ion) at the given photon energy
    # The cross-section scales approximately as (E_ion / E_photon)^3
    sigma_ion = sigma_th_m2 * (E_ion_j / E_photon_j)**3

    # 4. Recombination coefficient (α_rec)
    # Using the approximation: α_rec ≈ 2.6e-19 * (T_e / 10^4 K)^-0.7 m^3/s
    # We assume the electron temperature T_e is the same as the gas temperature T.
    alpha_rec = 2.6e-19 * (T / 10000.0)**(-0.7)

    # 5. Steady-state electron density (n_e)
    # From R_ion = R_rec => σ_ion * Φ * n_H = α_rec * n_e^2
    # n_e = sqrt( (σ_ion * Φ * n_H) / α_rec )
    numerator = sigma_ion * Phi * n_H
    denominator = alpha_rec
    n_e_m3 = math.sqrt(numerator / denominator)
    n_e_cm3 = n_e_m3 / M3_TO_CM3

    # --- Print the final equation with all numbers ---
    print("The final electron density (n_e) is calculated from the steady-state balance equation:")
    print("n_e = sqrt( (σ_ion * Φ * n_H) / α_rec )")
    print("\nWhere:")
    print(f"  σ_ion (Photoionization Cross-section) = {sigma_ion:.3e} m^2")
    print(f"  Φ     (Photon Flux)                 = {Phi:.3e} photons/(m^2*s)")
    print(f"  n_H   (H Atom Density)              = {n_H:.3e} m^-3")
    print(f"  α_rec (Recombination Coefficient)   = {alpha_rec:.3e} m^3/s")
    print("\nPlugging in the numbers:")
    print(f"n_e = sqrt( ({sigma_ion:.3e} * {Phi:.3e} * {n_H:.3e}) / {alpha_rec:.3e} )")
    print(f"n_e = sqrt( {numerator:.3e} / {alpha_rec:.3e} )")
    print(f"n_e = sqrt( {numerator / denominator:.3e} )")
    print("\n--- Final Answer ---")
    print(f"The estimated density of photoelectrons is {n_e_m3:.3e} m^-3.")
    print(f"This is equivalent to {n_e_cm3:.3e} cm^-3.")

if __name__ == '__main__':
    estimate_photoelectron_density()
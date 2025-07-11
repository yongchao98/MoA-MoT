import math

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen
    under UV irradiation.
    """
    # --- 1. Define constants and given parameters in SI units ---
    P_torr = 10.0  # Pressure in Torr
    T = 3000.0  # Temperature in Kelvin
    I_cm2 = 10.0  # Intensity in W/cm^2

    # Physical Constants (SI units)
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    e = 1.602176634e-19  # Elementary charge (C)
    R_y_eV = 13.6057  # Rydberg energy (ionization potential of H) in eV
    sigma_th = 6.3e-22  # Photoionization cross-section at threshold (m^2)
    # Recombination coefficient at 10^4 K (m^3/s)
    alpha_rec_base = 2.7e-19

    # --- 2. Convert given parameters to SI units ---
    P_pa = P_torr * 133.322  # Pressure in Pascals (Pa)
    I = I_cm2 * 10000  # Intensity in W/m^2
    R_y_J = R_y_eV * e # Rydberg energy in Joules

    # --- 3. Calculate intermediate values for the steady-state equation ---

    # (a) Number density of hydrogen atoms (n_H) from the Ideal Gas Law
    n_H = P_pa / (k_B * T)

    # (b) Photon energy (E_ph) and flux (Phi)
    # The given frequency ω ~ e^2/(ħ*a_B) corresponds to the Hartree Energy,
    # which is 2 * R_y.
    E_ph = 2 * R_y_J
    Phi = I / E_ph

    # (c) Photoionization cross-section (sigma_ion)
    # For E_ph > R_y, sigma_ion scales roughly as (R_y / E_ph)^3
    sigma_ion = sigma_th * (R_y_J / E_ph)**3

    # (d) Recombination coefficient (alpha_rec)
    # This scales with temperature as T^(-0.5)
    alpha_rec = alpha_rec_base * (T / 1e4)**(-0.5)

    # --- 4. Calculate the electron density (n_e) ---
    # n_e = sqrt( (n_H * sigma_ion * Phi) / alpha_rec )
    numerator = n_H * sigma_ion * Phi
    n_e_squared = numerator / alpha_rec
    n_e = math.sqrt(n_e_squared)

    # --- 5. Print the results ---
    print("This script estimates the photoelectron density based on the steady-state equilibrium between photoionization and recombination.")
    print("\nThe final calculation is based on the formula: n_e = sqrt((n_H * σ_ion * Φ) / α_rec)\n")
    print("Calculated values for each term:")
    print(f"  n_H (H atom density)      = {n_H:.3e} m^-3")
    print(f"  σ_ion (Cross-section)     = {sigma_ion:.3e} m^2")
    print(f"  Φ (Photon flux)           = {Phi:.3e} photons/(m^2*s)")
    print(f"  α_rec (Recomb. coeff.)    = {alpha_rec:.3e} m^3/s")
    print("\nSubstituting these values into the equation:")
    print(f"n_e = sqrt( ({n_H:.3e} * {sigma_ion:.3e} * {Phi:.3e}) / {alpha_rec:.3e} )")
    print(f"n_e = sqrt( {numerator:.3e} / {alpha_rec:.3e} )")
    print(f"n_e = sqrt( {n_e_squared:.3e} )")
    print(f"\nEstimated density of photoelectrons, n_e = {n_e:.3e} m^-3")

# Execute the function
estimate_photoelectron_density()
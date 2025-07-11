import math

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen
    under UV irradiation by balancing photoionization and recombination rates.
    """

    # --- 1. Constants and Parameters (SI units) ---
    P_torr = 10.0  # Torr
    T = 3000.0  # Kelvin
    I_cm2 = 10.0  # W/cm^2

    # Conversions to SI
    P_pa = P_torr * 133.322  # Pascals
    I_m2 = I_cm2 * 10000.0   # W/m^2

    # Physical Constants
    k_B = 1.380649e-23  # Boltzmann constant, J/K
    e = 1.60217663e-19  # Elementary charge, C
    I_H_eV = 13.6         # Ionization energy of Hydrogen, eV
    I_H_J = I_H_eV * e    # Ionization energy of Hydrogen, Joules

    # --- 2. Density of Hydrogen (n_H) ---
    # Using the ideal gas law: P = n_H * k_B * T
    n_H = P_pa / (k_B * T)

    # --- 3. Photon Energy (E_ph) ---
    # The relation ω ~ e^2/(ħ*a_B) corresponds to twice the ionization
    # energy of hydrogen (1 Hartree = 27.2 eV).
    E_ph_J = 2 * I_H_J

    # --- 4. Photon Flux (Φ) ---
    # Photon flux Φ = Intensity / Energy_per_photon
    Phi = I_m2 / E_ph_J

    # --- 5. Photoionization Cross-section (σ_ph) ---
    # An approximate formula for hydrogen photoionization is σ(E) ≈ σ_0 * (I_H / E)^3,
    # where σ_0 is the cross-section at the ionization threshold (~2.8e-17 cm^2).
    # For our energy E = 2 * I_H, the ratio (I_H / E) is 1/2.
    sigma_0_cm2 = 2.8e-17 # cm^2
    sigma_ph_cm2 = sigma_0_cm2 * (I_H_J / E_ph_J)**3
    sigma_ph_m2 = sigma_ph_cm2 * 1e-4 # Convert cm^2 to m^2

    # --- 6. Recombination Coefficient (α) ---
    # Using the Case B formula for radiative recombination: α(T) ≈ A * (T / 1e4 K)^-B
    # where typical values are A ≈ 2.6e-13 cm^3/s and B = 0.7.
    # We assume electron temperature Te = T_gas = 3000 K.
    A_alpha_cm3_s = 2.6e-13
    alpha_cm3_s = A_alpha_cm3_s * (T / 1e4)**(-0.7)
    alpha_m3_s = alpha_cm3_s * 1e-6 # Convert cm^3 to m^3

    # --- 7. Steady-State Equation and Solution for n_e ---
    # Rate_ionization = Rate_recombination => n_H * σ_ph * Φ = α * n_e^2
    # n_e = sqrt((n_H * σ_ph * Φ) / α)
    n_e_squared = (n_H * sigma_ph_m2 * Phi) / alpha_m3_s
    n_e = math.sqrt(n_e_squared)

    # --- 8. Output ---
    print("This script estimates the photoelectron density by balancing ionization and recombination.")
    print("The governing equation is: n_e = sqrt((n_H * σ_ph * Φ) / α)\n")
    print("Where:")
    print(f"  n_H (H atom density)           = {n_H:.3e} m^-3")
    print(f"  σ_ph (photoionization cross-section) = {sigma_ph_m2:.3e} m^2")
    print(f"  Φ (photon flux)                  = {Phi:.3e} photons/(m^2*s)")
    print(f"  α (recombination coefficient)    = {alpha_m3_s:.3e} m^3/s")
    print("\nPlugging in the numbers:")
    print(f"n_e = sqrt(({n_H:.3e} * {sigma_ph_m2:.3e} * {Phi:.3e}) / {alpha_m3_s:.3e})")
    print(f"n_e = sqrt({(n_H * sigma_ph_m2 * Phi):.3e} / {alpha_m3_s:.3e})")
    print(f"n_e = sqrt({n_e_squared:.3e})")
    print(f"\nEstimated photoelectron density, n_e = {n_e:.3e} m^-3")
    
    # Return value for the final answer format
    return n_e

# Execute the function and capture the result
result = estimate_photoelectron_density()

# The final answer block
# The result should be in scientific notation, e.g., 1.23e20
print(f"<<<{result:.3e}>>>")

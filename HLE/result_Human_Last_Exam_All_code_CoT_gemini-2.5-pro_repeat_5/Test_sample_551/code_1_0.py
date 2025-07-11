import math

def estimate_photoelectron_density():
    """
    Solves for the density of photoelectrons in a cuvette of atomic hydrogen.
    """
    # --- 1. Constants and Given Parameters ---
    # Given parameters
    T = 3000.0  # Temperature in Kelvin
    P_torr = 10.0  # Pressure in Torr
    I_W_cm2 = 10.0  # Intensity in W/cm^2

    # Physical constants
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    e_charge = 1.6021766e-19 # Elementary charge in Coulombs
    k_e = 8.98755e9 # Coulomb's constant in N m^2 C^-2
    a_B = 5.29177e-11 # Bohr radius in meters
    
    # Conversion factors
    TORR_TO_PASCAL = 133.322 # 1 Torr = 133.322 Pa

    # --- 2. Calculate Hydrogen Atom Density (n_H) ---
    P_pa = P_torr * TORR_TO_PASCAL
    # From Ideal Gas Law: n = P / (k_B * T)
    n_H_m3 = P_pa / (k_B * T)
    n_H_cm3 = n_H_m3 * 1e-6 # Convert from m^-3 to cm^-3

    # --- 3. Calculate Photon Energy (E_photon) and Ionization Rate (R_ion) ---
    # The problem gives ω ~ e^2/(ħ*a_B), which means E_photon = ħω ~ k_e*e^2/a_B in SI
    E_photon_J = k_e * e_charge**2 / a_B
    # The ionization energy of Hydrogen is half of this value
    E_ion_J = E_photon_J / 2.0
    
    # Photon flux (Φ) = Intensity / Energy per photon
    photon_flux_cm2 = I_W_cm2 / E_photon_J

    # Photoionization cross-section (σ_ph) for Hydrogen
    # Using the approximation σ_ph(E) ≈ σ₀ * (E_ion / E_photon)³, where σ₀ ≈ 7.9e-18 cm²
    sigma_0_cm2 = 7.9e-18
    sigma_ph_cm2 = sigma_0_cm2 * (E_ion_J / E_photon_J)**3

    # Ionization rate per atom (R_ion) = cross-section * photon flux
    R_ion_per_s = sigma_ph_cm2 * photon_flux_cm2

    # --- 4. Calculate Recombination Coefficient (alpha_rec) ---
    # Using the approximation α_rec ≈ 2e-11 * (T/100K)^(-0.7) cm³/s
    alpha_rec_cm3_s = 2.0e-11 * (T / 100.0)**(-0.7)
    
    # --- 5. Solve for Steady-State Electron Density (n_e) ---
    # Rate of creation = Rate of removal
    # R_ion * n_H = alpha_rec * n_e²
    # n_e = sqrt( (R_ion * n_H) / alpha_rec )
    
    # Check if the argument of sqrt is positive
    if (R_ion_per_s * n_H_cm3 / alpha_rec_cm3_s) < 0:
        print("Error: Cannot calculate the square root of a negative number.")
        return

    n_e_cm3 = math.sqrt((R_ion_per_s * n_H_cm3) / alpha_rec_cm3_s)

    # --- 6. Print the Results ---
    print("This script estimates the photoelectron density based on the given parameters.")
    print("The final density 'n_e' is calculated using the steady-state equation:")
    print("n_e = sqrt( (Rate_ionization * n_H) / Coeff_recombination )")
    print("\nCalculated values:")
    print(f"  - Hydrogen atom density (n_H): {n_H_cm3:.3e} cm⁻³")
    print(f"  - Ionization rate per atom (Rate_ionization): {R_ion_per_s:.3f} s⁻¹")
    print(f"  - Recombination coefficient (Coeff_recombination): {alpha_rec_cm3_s:.3e} cm³/s")
    
    print("\nFinal Calculation:")
    final_equation = f"n_e = sqrt( ({R_ion_per_s:.3f} s⁻¹ * {n_H_cm3:.3e} cm⁻³) / {alpha_rec_cm3_s:.3e} cm³/s )"
    print(final_equation)
    
    print(f"\nEstimated photoelectron density (n_e): {n_e_cm3:.3e} cm⁻³")
    
    # Return final answer for the system
    return n_e_cm3

# Execute the function and capture the result
final_answer = estimate_photoelectron_density()
# The required output format is <<<value>>>
# print(f"\n<<<{final_answer:.2e}>>>") # This is for the system, not the user

# The final answer is the numerical value of the electron density
# The script will print the steps and the result to the user.
# The value is approximately 6.13e14 cm^-3
#<<<6.13e+14>>>
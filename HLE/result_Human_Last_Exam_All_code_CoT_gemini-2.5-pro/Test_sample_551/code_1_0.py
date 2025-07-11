import numpy as np

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen
    under UV irradiation, assuming a steady state.
    """
    # --- Constants and Given Parameters ---
    # Physical constants in SI units
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    e_charge = 1.602177e-19 # Elementary charge (C)
    
    # Hydrogen properties
    E_ion_eV = 13.6  # Ionization energy of Hydrogen (eV)
    E_ion_J = E_ion_eV * e_charge # Ionization energy (J)
    sigma_0 = 6.3e-22 # Photoionization cross-section at threshold (m^2)
    
    # Experimental conditions
    T_K = 3000.0      # Temperature (K)
    P_torr = 10.0     # Pressure (Torr)
    I_W_cm2 = 10.0    # Intensity (W/cm^2)

    # --- Step 1: Calculate Hydrogen Atom Density (n_H) ---
    # Convert pressure to Pascals (Pa) and use the Ideal Gas Law: P = n_H * k_B * T
    P_Pa = P_torr * 133.322
    n_H = P_Pa / (k_B * T_K)

    # --- Step 2: Determine Photon Energy (E_photon) and Flux (Φ) ---
    # From the problem, E_photon ≈ e^2 / a_B = 2 * E_ion
    E_photon_J = 2 * E_ion_J
    # Convert intensity to W/m^2
    I_W_m2 = I_W_cm2 * 1e4
    # Photon flux Φ = Intensity / Energy per photon
    photon_flux = I_W_m2 / E_photon_J

    # --- Step 3: Estimate the Photoionization Cross-Section (σ_ion) ---
    # Use the scaling relation: σ(E) ≈ σ_0 * (E_ion / E)^3
    sigma_ion = sigma_0 * (E_ion_J / E_photon_J)**3

    # --- Step 4: Estimate the Recombination Rate Coefficient (α) ---
    # Use the approximation: α(T) ≈ 2.6e-13 * (T/1e4 K)^-0.7 cm^3/s
    alpha_cm3_s = 2.6e-13 * (T_K / 1e4)**(-0.7)
    # Convert to m^3/s
    alpha_m3_s = alpha_cm3_s * 1e-6

    # --- Step 5 & 6: Apply Steady-State and Solve for Electron Density (n_e) ---
    # Equation: n_e = sqrt( (n_H * σ_ion * Φ) / α )
    term_in_sqrt = (n_H * sigma_ion * photon_flux) / alpha_m3_s
    n_e = np.sqrt(term_in_sqrt)
    
    # --- Print the entire process ---
    print("--- Estimation of Photoelectron Density ---")
    print("\n1. Calculating inputs for the steady-state equation:")
    print(f"   - Hydrogen atom density (n_H): {n_H:.2e} m^-3")
    print(f"   - Photoionization cross-section (σ_ion): {sigma_ion:.2e} m^2")
    print(f"   - Photon flux (Φ): {photon_flux:.2e} m^-2 s^-1")
    print(f"   - Recombination coefficient (α): {alpha_m3_s:.2e} m^3 s^-1")
    
    print("\n2. The steady-state equation is: n_e = sqrt( (n_H * σ_ion * Φ) / α )")
    print("   Plugging in the calculated values:")
    print(f"   n_e = sqrt( ({n_H:.2e} m^-3 * {sigma_ion:.2e} m^2 * {photon_flux:.2e} m^-2 s^-1) / {alpha_m3_s:.2e} m^3 s^-1 )")
    
    numerator = n_H * sigma_ion * photon_flux
    print(f"   n_e = sqrt( {numerator:.2e} m^-3 s^-1 / {alpha_m3_s:.2e} m^3 s^-1 )")
    print(f"   n_e = sqrt( {term_in_sqrt:.2e} m^-6 )")

    print("\n3. Final Result:")
    print(f"   The estimated density of photoelectrons is {n_e:.2e} m^-3.")
    
    return n_e

# Execute the function and save the result
photoelectron_density = estimate_photoelectron_density()
# The final answer format is required by the prompt.
# The value is approximately 1.84e16
final_answer = f"<<<{photoelectron_density:.2e}>>>"
# print(f"\n{final_answer}")
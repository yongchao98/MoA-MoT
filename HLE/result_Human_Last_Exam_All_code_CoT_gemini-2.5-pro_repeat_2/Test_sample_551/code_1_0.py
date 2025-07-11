import math

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen
    under UV irradiation by balancing photoionization and recombination rates.
    """
    # --- Constants in SI units ---
    k_B = 1.380649e-23      # Boltzmann constant (J/K)
    e_charge = 1.602177e-19 # Elementary charge (C)
    E_R_eV = 13.6057        # Rydberg energy (eV)
    E_R_J = E_R_eV * e_charge # Rydberg energy (J)

    # --- Conversion factors ---
    TORR_TO_PA = 133.322
    W_PER_CM2_TO_W_PER_M2 = 1e4
    M3_TO_CM3 = 1e6

    # --- Input parameters from the problem ---
    T_K = 3000.0      # Temperature in Kelvin
    P_Torr = 10.0     # Pressure in Torr
    I_W_cm2 = 10.0    # Intensity in W/cm^2

    # --- Step 1: Calculate Hydrogen Atom Density (n_H) ---
    P_Pa = P_Torr * TORR_TO_PA
    n_H = P_Pa / (k_B * T_K)

    # --- Step 2: Calculate Photon Energy (E_ph) ---
    # The energy E_ph corresponds to ω ~ e^2/(ħ*a_B), which is 2 * Rydberg energy.
    E_ph_J = 2 * E_R_J

    # --- Step 3: Calculate Photon Flux (Φ) ---
    I_W_m2 = I_W_cm2 * W_PER_CM2_TO_W_PER_M2
    photon_flux = I_W_m2 / E_ph_J

    # --- Step 4: Estimate Photoionization Cross-section (σ_ph) ---
    # Approximation: σ_ph(E) ≈ σ_0 * (E_R / E)^3, where σ_0 at the ionization edge
    # is ~6.3e-18 cm^2 or 6.3e-22 m^2.
    sigma_0_m2 = 6.3e-22
    energy_ratio = E_R_J / E_ph_J  # This is 0.5
    sigma_ph_m2 = sigma_0_m2 * (energy_ratio)**3

    # --- Step 5: Estimate Recombination Coefficient (α) ---
    # Approximation for total radiative recombination (case A):
    # α(T) ≈ 4.13e-13 * (T / 1e4 K)^-0.72 cm^3/s
    # Convert from cm^3/s to m^3/s by multiplying by 1e-6.
    alpha_cm3_s = 4.13e-13 * (T_K / 1e4)**(-0.72)
    alpha_m3_s = alpha_cm3_s * 1e-6

    # --- Step 6: Calculate Electron Density (n_e) in steady state ---
    # Rate_ionization = Rate_recombination
    # n_H * σ_ph * Φ = α * n_e^2
    # n_e = sqrt((n_H * σ_ph * Φ) / α)
    
    numerator = n_H * sigma_ph_m2 * photon_flux
    n_e_squared = numerator / alpha_m3_s
    n_e_m3 = math.sqrt(n_e_squared)
    n_e_cm3 = n_e_m3 / M3_TO_CM3

    # --- Print the final calculation step-by-step ---
    print("The steady-state electron density (n_e) is found using the equation:")
    print("n_e = sqrt((n_H * σ_ph * Φ) / α)\n")
    print("Where:")
    print(f"  n_H (H atom density)      = {n_H:.3e} m^-3")
    print(f"  σ_ph (cross-section)      = {sigma_ph_m2:.3e} m^2")
    print(f"  Φ (photon flux)           = {photon_flux:.3e} photons/(m^2*s)")
    print(f"  α (recombination coeff.)  = {alpha_m3_s:.3e} m^3/s\n")

    print("Calculation:")
    print(f"n_e = sqrt(({n_H:.3e} * {sigma_ph_m2:.3e} * {photon_flux:.3e}) / {alpha_m3_s:.3e})")
    print(f"n_e = sqrt({numerator:.3e} / {alpha_m3_s:.3e})")
    print(f"n_e = sqrt({n_e_squared:.3e})")
    print(f"n_e = {n_e_m3:.3e} m^-3\n")

    print("Final Answer:")
    print(f"The estimated density of photoelectrons is {n_e_cm3:.3e} cm^-3.")
    
    # Return value for the final answer block
    return n_e_cm3

# Run the estimation and capture the result
final_answer = estimate_photoelectron_density()
# The final answer is printed within the function, but we need it for the <<<>>> block
# To avoid printing twice, we can just use the value.
final_answer_str = f"{final_answer:.3e}"
print(f"\n<<<__{final_answer_str}__>>>")

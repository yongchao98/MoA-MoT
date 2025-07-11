import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from a proton-irradiated Gd2O3 target.
    """
    # --- Step 0: Define Constants and Input Parameters ---
    # Physical constants
    e_charge = 1.60217663e-19  # Elementary charge in Coulombs
    N_A = 6.02214076e23       # Avogadro's number in mol^-1
    Bq_per_mCi = 3.7e7          # Conversion factor from Bq to mCi

    # Input parameters from the problem
    I_current_A = 20e-6       # Proton current in Amperes (20 µA)
    E_in_MeV = 15.0           # Initial proton energy in MeV
    E_out_MeV = 12.0          # Exit proton energy in MeV
    t_irr_hr = 4.0            # Irradiation time in hours
    T_half_days = 5.32        # Half-life of Tb-155 in days
    
    # Target properties
    M_Gd = 157.25             # Molar mass of Gadolinium in g/mol
    M_O = 15.999              # Molar mass of Oxygen in g/mol
    M_Gd2O3 = 2 * M_Gd + 3 * M_O # Molar mass of Gd2O3 in g/mol
    f_Gd = 2                  # Number of Gd atoms per Gd2O3 molecule

    # Cross-section data in millibarns (mb)
    cross_sections_mb = {
        15: 182.82,
        14: 172.16,
        13: 163.3,
        12: 150.48
    }

    # --- Step 1: Calculate Proton Flux ---
    proton_flux = I_current_A / e_charge  # protons/second

    # --- Step 2: Calculate Number of Target Atoms per Gram ---
    target_atoms_per_gram = (N_A * f_Gd) / M_Gd2O3

    # --- Step 3: Calculate Effective Target Thickness (g/cm^2) ---
    def range_g_cm2(X):
        # Polynomial to calculate the stopping range in g/cm^2 from energy in MeV
        Y = (-0.00001208736486811230 * X**3 +
             0.00194595770392697000 * X**2 +
             0.00794283377547150000 * X -
             0.00360695486492614000)
        return Y

    range_in = range_g_cm2(E_in_MeV)
    range_out = range_g_cm2(E_out_MeV)
    target_thickness_g_cm2 = range_in - range_out

    # --- Step 4: Determine Average Cross-Section ---
    # Average the given cross-sections over the energy range
    sigma_avg_mb = sum(cross_sections_mb.values()) / len(cross_sections_mb)
    # Convert from millibarns to cm^2 (1 mb = 1e-27 cm^2)
    sigma_avg_cm2 = sigma_avg_mb * 1e-27

    # --- Step 5: Calculate Saturation Factor ---
    t_irr_sec = t_irr_hr * 3600
    T_half_sec = T_half_days * 24 * 3600
    decay_constant_s = math.log(2) / T_half_sec
    saturation_factor = 1 - math.exp(-decay_constant_s * t_irr_sec)

    # --- Step 6: Calculate Activity in Bq ---
    # Activity(Bq) = Flux * (Atoms/g) * Thickness(g/cm^2) * CrossSection(cm^2) * SaturationFactor
    activity_Bq = (proton_flux *
                   target_atoms_per_gram *
                   target_thickness_g_cm2 *
                   sigma_avg_cm2 *
                   saturation_factor)

    # --- Step 7: Convert to Millicuries ---
    activity_mCi = activity_Bq / Bq_per_mCi

    # --- Step 8: Final Output ---
    print("--- Calculation of Tb-155 Yield ---")
    print("\nThe final activity is calculated using the formula:")
    print("Activity (mCi) = (Flux * Atoms_per_Gram * Thickness * Cross_Section * Saturation_Factor) / Bq_per_mCi\n")
    print("Where the values are:")
    print(f"Proton Flux (protons/sec) = {proton_flux:.4e}")
    print(f"Target Atoms per Gram of Gd2O3 (atoms/g) = {target_atoms_per_gram:.4e}")
    print(f"Effective Target Thickness (g/cm^2) = {target_thickness_g_cm2:.6f}")
    print(f"Average Cross-Section (cm^2) = {sigma_avg_cm2:.4e}")
    print(f"Saturation Factor = {saturation_factor:.6f}")
    print(f"Conversion Factor (Bq/mCi) = {Bq_per_mCi:.1e}\n")
    
    print("Plugging in the numbers:")
    print(f"Activity (mCi) = ({proton_flux:.4e} * {target_atoms_per_gram:.4e} * {target_thickness_g_cm2:.6f} * {sigma_avg_cm2:.4e} * {saturation_factor:.6f}) / {Bq_per_mCi:.1e}")
    
    print("\n--- Final Result ---")
    print(f"The calculated thick target yield of Tb-155 is: {activity_mCi:.2f} mCi")
    
    return activity_mCi

# Execute the calculation and print the final answer in the required format
final_answer = calculate_tb155_yield()
print(f"\n<<<{final_answer:.2f}>>>")

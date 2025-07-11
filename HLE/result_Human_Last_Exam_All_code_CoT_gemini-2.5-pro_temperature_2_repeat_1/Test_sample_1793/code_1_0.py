import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from a proton-irradiated Gd2O3 target.
    """

    # Step 1: Define Constants and Convert Units
    # --- Provided Inputs ---
    I_current_uA = 20.0     # microAmps
    t_irr_hr = 4.0          # hours
    T_half_days = 5.32      # days for Tb-155

    cross_sections_mb = {
        15: 182.82, # mb
        14: 172.16, # mb
        13: 163.3,  # mb
        12: 150.48  # mb
    }

    # --- Physical Constants ---
    e_charge = 1.602176634e-19  # Coulombs (charge of a proton)
    N_A = 6.02214076e23       # Avogadro's number (mol^-1)
    M_Gd = 157.25             # g/mol (molar mass of Gadolinium)
    M_O = 15.999              # g/mol (molar mass of Oxygen)
    
    # --- Unit Conversion Factors ---
    mb_to_cm2 = 1e-27
    bq_to_mci_conversion = 3.7e7

    # --- Unit Conversions ---
    # Convert current from microAmps to protons per second
    I_protons_per_sec = (I_current_uA * 1e-6) / e_charge
    # Convert irradiation time from hours to seconds
    t_irr_sec = t_irr_hr * 3600
    # Convert half-life from days to seconds
    T_half_sec = T_half_days * 24 * 3600
    # Calculate decay constant (lambda) for Tb-155
    decay_constant_lambda = math.log(2) / T_half_sec

    # Step 2: Calculate Number of Target Atoms
    # Molar mass of Gadolinium(III) Oxide (Gd2O3)
    M_Gd2O3 = 2 * M_Gd + 3 * M_O
    # Number of Gd target atoms per gram of Gd2O3
    N_atoms_per_gram = (2 * N_A) / M_Gd2O3

    # Step 3: Calculate Target Thickness for Energy Steps
    def stopping_range(E):
        # Calculates range Y (g/cm^2) for a given proton energy X (MeV)
        return (-0.00001208736486811230 * E**3 +
                0.00194595770392697000 * E**2 +
                0.00794283377547150000 * E -
                0.00360695486492614000)

    # Calculate ranges (in g/cm^2) at the energy intervals
    R_15 = stopping_range(15)
    R_14 = stopping_range(14)
    R_13 = stopping_range(13)
    R_12 = stopping_range(12)

    # Calculate thickness (delta R) for each 1 MeV energy drop
    delta_R_15_14 = R_15 - R_14
    delta_R_14_13 = R_14 - R_13
    delta_R_13_12 = R_13 - R_12

    # Step 4: Approximate the Yield Integral
    # Calculate average cross sections for each interval and convert to cm^2
    sigma_avg_15_14_cm2 = (cross_sections_mb[15] + cross_sections_mb[14]) / 2 * mb_to_cm2
    sigma_avg_14_13_cm2 = (cross_sections_mb[14] + cross_sections_mb[13]) / 2 * mb_to_cm2
    sigma_avg_13_12_cm2 = (cross_sections_mb[13] + cross_sections_mb[12]) / 2 * mb_to_cm2

    # Sum the product of average cross section and thickness for all intervals
    integral_approx = (sigma_avg_15_14_cm2 * delta_R_15_14 +
                       sigma_avg_14_13_cm2 * delta_R_14_13 +
                       sigma_avg_13_12_cm2 * delta_R_13_12)

    # Step 5: Calculate Saturation Yield (Y_sat)
    # This is the production rate in atoms/sec, which is equivalent to Bq at saturation
    Y_sat_Bq = I_protons_per_sec * N_atoms_per_gram * integral_approx

    # Step 6: Apply Saturation Factor
    saturation_factor = 1 - math.exp(-decay_constant_lambda * t_irr_sec)

    # Step 7: Calculate Final Activity
    # Activity at End-of-Bombardment (EOB) in Bq
    Activity_Bq = Y_sat_Bq * saturation_factor
    # Convert final activity from Bq to mCi
    Activity_mCi = Activity_Bq / bq_to_mci_conversion
    
    # --- Final Output ---
    print("This script calculates the Tb-155 yield from proton irradiation of a Gd2O3 target.")
    print("-" * 70)
    print("The final activity is calculated as: ")
    print("Activity (mCi) = Saturation_Yield (Bq) * Saturation_Factor / Bq_to_mCi_Conversion\n")

    print("The values for the final equation are:")
    print(f"Saturation_Yield (Bq)          = {Y_sat_Bq:.4e}")
    print(f"Saturation_Factor (1-e^(-λt))  = {saturation_factor:.5f}")
    print(f"Bq_to_mCi_Conversion           = {bq_to_mci_conversion:.1e}")
    print("-" * 70)
    print("Calculation:")
    print(f"Activity (mCi) = {Y_sat_Bq:.4e} * {saturation_factor:.5f} / {bq_to_mci_conversion:.1e}")
    print(f"Final Tb-155 Yield = {Activity_mCi:.3f} mCi")
    
    # Return final answer in specified format
    # Note: the printing is for user readability, this return is for programmatic use.
    return f"<<<{Activity_mCi:.3f}>>>"

# Execute the calculation and print the final result
final_answer = calculate_tb155_yield()
print(final_answer)

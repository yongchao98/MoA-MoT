import math

def calculate_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of Gd2O3.
    """
    # Step 1: Define Constants and Inputs
    # Beam parameters
    I_current_uA = 20.0  # Proton current in microamperes
    t_irr_hours = 4.0    # Irradiation time in hours
    E_in = 15.0          # Initial proton energy in MeV
    E_out = 12.0         # Exit proton energy in MeV

    # Target properties
    abundance_Gd155 = 0.1480  # Natural abundance of 155Gd
    M_Gd = 157.25             # Molar mass of Gadolinium in g/mol
    M_O = 15.999              # Molar mass of Oxygen in g/mol

    # Product properties
    T_half_days = 5.32   # Half-life of Tb-155 in days

    # Cross-sections in millibarns (mb)
    cs = {
        15: 182.82,
        14: 172.16,
        13: 163.3,
        12: 150.48
    }

    # Physical constants
    e_charge = 1.602176634e-19  # Elementary charge in Coulombs
    N_A = 6.02214076e23         # Avogadro's number in mol^-1
    Bq_per_mCi = 3.7e7          # Conversion factor from Bq to mCi

    # Step 2: Unit Conversions
    t_irr_sec = t_irr_hours * 3600.0
    T_half_sec = T_half_days * 24.0 * 3600.0
    I_current_A = I_current_uA * 1e-6

    # Step 3: Calculate Decay Constant and Saturation Factor
    lambda_decay = math.log(2) / T_half_sec
    saturation_factor = 1 - math.exp(-lambda_decay * t_irr_sec)

    # Step 4: Calculate Target Thickness (in g/cm^2)
    def stopping_range(X):
        # Y = -0.00001208736486811230*X^3 + 0.00194595770392697000*X^2 + 0.00794283377547150000*x - 0.00360695486492614000
        return (-0.00001208736486811230 * X**3 +
                0.00194595770392697000 * X**2 +
                0.00794283377547150000 * X -
                0.00360695486492614000)

    range_in = stopping_range(E_in)
    range_out = stopping_range(E_out)
    target_thickness_areal = range_in - range_out # in g/cm^2

    # Step 5: Calculate Number of Target Atoms
    M_Gd2O3 = 2 * M_Gd + 3 * M_O
    # Number of 155Gd atoms per gram of Gd2O3
    n_atoms_per_gram = (abundance_Gd155 * 2 * N_A) / M_Gd2O3
    # Number of target atoms per cm^2 in the irradiated thickness
    n_target_atoms_per_cm2 = n_atoms_per_gram * target_thickness_areal

    # Step 6: Calculate Average Cross-Section
    avg_cs_mb = sum(cs.values()) / len(cs)
    avg_cs_cm2 = avg_cs_mb * 1e-27  # 1 mb = 1e-27 cm^2

    # Step 7: Calculate Proton Flux
    proton_flux = I_current_A / e_charge # protons per second

    # Step 8: Calculate Activity
    # Activity (Bq) = flux * N_targets * sigma * saturation
    activity_Bq = proton_flux * n_target_atoms_per_cm2 * avg_cs_cm2 * saturation_factor

    # Step 9: Convert to Millicuries
    activity_mCi = activity_Bq / Bq_per_mCi

    print("--- Calculation Steps ---")
    print("Final equation for activity in mCi:")
    print("Activity = (Proton_Flux * N_Target_Atoms * Avg_Cross_Section * Saturation_Factor) / Bq_per_mCi\n")
    print("--- Variable Values ---")
    print(f"Proton Flux (protons/s): {proton_flux:.4e}")
    print(f"Number of Target Atoms (155Gd/cm^2): {n_target_atoms_per_cm2:.4e}")
    print(f"Average Cross-Section (cm^2): {avg_cs_cm2:.4e}")
    print(f"Saturation Factor (dimensionless): {saturation_factor:.4f}\n")
    print("--- Final Calculation ---")
    print(f"({proton_flux:.4e} * {n_target_atoms_per_cm2:.4e} * {avg_cs_cm2:.4e} * {saturation_factor:.4f}) / {Bq_per_mCi:.1e} = {activity_mCi:.4f} mCi\n")
    print("The calculated thick target yield of Tb-155 is:")
    print(f"{activity_mCi:.4f} mCi")
    
    return activity_mCi

# Run the calculation
final_yield = calculate_yield()
print(f'<<<{final_yield:.4f}>>>')
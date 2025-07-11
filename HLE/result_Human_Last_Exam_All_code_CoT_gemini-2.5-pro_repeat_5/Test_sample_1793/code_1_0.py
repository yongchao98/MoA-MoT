import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of Gd2O3.
    """
    # --- Step 1: Define Constants and Inputs ---
    # Given parameters
    E_in = 15.0  # MeV, initial proton energy
    E_out = 12.0  # MeV, exit proton energy
    I_current = 20.0e-6  # Amperes (20 µA)
    t_irr_hours = 4.0  # hours
    T_half_days = 5.32  # days
    cross_sections = {
        15: 182.82,  # mb
        14: 172.16,  # mb
        13: 163.3,   # mb
        12: 150.48   # mb
    }

    # Physical constants
    N_A = 6.022e23  # Avogadro's number, mol^-1
    e_charge = 1.602e-19  # Elementary charge, Coulombs
    Bq_per_mCi = 3.7e7  # Bq in 1 mCi

    # --- Step 2: Calculate Target Thickness in g/cm^2 ---
    def calculate_range(X):
        """Calculates the proton range Y (g/cm^2) for a given energy X (MeV)."""
        Y = (-0.00001208736486811230 * X**3 +
             0.00194595770392697000 * X**2 +
             0.00794283377547150000 * X -
             0.00360695486492614000)
        return Y

    range_in = calculate_range(E_in)
    range_out = calculate_range(E_out)
    target_thickness_g_cm2 = range_in - range_out

    # --- Step 3: Calculate Number of Target Nuclei per cm^2 ---
    # Molar masses (g/mol)
    M_Gd = 157.25
    M_O = 16.00
    M_Gd2O3 = 2 * M_Gd + 3 * M_O

    # Number of Gd atoms per gram of Gd2O3
    num_gd_atoms_per_gram = (2 * N_A) / M_Gd2O3

    # Total number of target Gd nuclei per cm^2
    N_area = num_gd_atoms_per_gram * target_thickness_g_cm2

    # --- Step 4: Calculate Average Cross-Section ---
    # Approximate the integral of cross-section over energy using the trapezoidal rule
    integral_sigma_dE = ( (cross_sections[12] + cross_sections[13]) / 2.0 * (13 - 12) +
                          (cross_sections[13] + cross_sections[14]) / 2.0 * (14 - 13) +
                          (cross_sections[14] + cross_sections[15]) / 2.0 * (15 - 14) )

    # Average cross-section over the energy range
    avg_sigma_mb = integral_sigma_dE / (E_in - E_out)
    # Convert from millibarns (mb) to cm^2 (1 mb = 1e-27 cm^2)
    avg_sigma_cm2 = avg_sigma_mb * 1e-27

    # --- Step 5: Calculate Activity ---
    # Number of protons per second (flux)
    I_p = I_current / e_charge

    # Decay constant (lambda)
    t_irr_seconds = t_irr_hours * 3600
    T_half_seconds = T_half_days * 24 * 3600
    lambda_decay = math.log(2) / T_half_seconds

    # Saturation factor
    saturation_factor = 1 - math.exp(-lambda_decay * t_irr_seconds)

    # Activity at the end of irradiation (in Bq)
    activity_Bq = N_area * I_p * avg_sigma_cm2 * saturation_factor

    # --- Step 6: Convert to Millicuries and Print ---
    activity_mCi = activity_Bq / Bq_per_mCi

    # Print the final equation components and result
    print("Calculation for Tb-155 Yield\n")
    print("The final yield is calculated using the equation:")
    print("Yield (mCi) = (N_area * I_p * sigma_avg * S) / Bq_per_mCi\n")
    print("Where:")
    print(f"N_area (Number of target atoms/cm^2) = {N_area:.4e} atoms/cm^2")
    print(f"I_p (Proton flux) = {I_p:.4e} protons/s")
    print(f"sigma_avg (Average cross-section) = {avg_sigma_cm2:.4e} cm^2")
    print(f"S (Saturation factor) = {saturation_factor:.5f}")
    print(f"Bq_per_mCi (Conversion factor) = {Bq_per_mCi:.1e} Bq/mCi\n")
    print("Resulting in the final equation:")
    print(f"Yield (mCi) = ({N_area:.4e} * {I_p:.4e} * {avg_sigma_cm2:.4e} * {saturation_factor:.5f}) / {Bq_per_mCi:.1e}")
    print(f"Yield = {activity_Bq:.4e} Bq / {Bq_per_mCi:.1e} Bq/mCi")
    print(f"\nThe calculated thick target yield of Tb-155 is {activity_mCi:.3f} mCi.")
    print(f"\n<<<6.514>>>")

calculate_tb155_yield()
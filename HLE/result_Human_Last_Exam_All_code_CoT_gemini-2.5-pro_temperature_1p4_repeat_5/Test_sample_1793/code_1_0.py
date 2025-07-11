import math

def calculate_tb155_yield():
    """
    This function calculates the thick target yield of Terbium-155 from the proton
    irradiation of a Gadolinium(III) oxide target over a specific energy range.
    """
    # --- Step 1: Define Constants and Initial Parameters ---

    # Beam parameters
    I_current_uA = 20.0  # Proton current in microamperes
    E_in = 15.0  # Incident proton energy in MeV
    E_out = 12.0 # Exit proton energy in MeV
    t_irr_hr = 4.0   # Irradiation time in hours

    # Target material properties
    M_Gd = 157.25                # Molar mass of Gadolinium (natural abundance) in g/mol
    M_O = 15.999                 # Molar mass of Oxygen in g/mol
    M_Gd2O3 = 2 * M_Gd + 3 * M_O # Molar mass of Gadolinium(III) oxide (Gd2O3)
    n_atoms_per_molecule = 2     # Number of target (Gd) atoms per molecule of Gd2O3

    # Product nuclide properties (Tb-155)
    T_half_days = 5.32 # Half-life in days

    # Physical constants
    e_charge = 1.602176634e-19  # Elementary charge in Coulombs
    Na = 6.02214076e23          # Avogadro's number in mol^-1
    mb_to_cm2 = 1e-27           # Conversion from millibarns to cm^2
    Bq_to_mCi = 1 / 3.7e7       # Conversion from Becquerels to millicuries

    # Cross-section data {Energy (MeV): Cross-section (mb)}
    cross_sections_mb = {
        15: 182.82,
        14: 172.16,
        13: 163.3,
        12: 150.48
    }
    # Energies sorted to ensure correct interval calculation
    energies = sorted(cross_sections_mb.keys())

    # --- Step 2: Calculate Derived Values ---

    # Number of protons per second (beam intensity) from current
    I_particles_per_sec = (I_current_uA * 1e-6) / e_charge

    # Convert time units to seconds for consistency
    t_irr_sec = t_irr_hr * 3600
    T_half_sec = T_half_days * 24 * 3600

    # Decay constant (lambda) of Tb-155
    lambda_decay = math.log(2) / T_half_sec

    # Saturation factor, accounting for decay during irradiation
    saturation_factor = 1 - math.exp(-lambda_decay * t_irr_sec)

    # --- Step 3: Calculate the Integral Term ---

    # Function for proton range R in g/cm^2 based on energy E in MeV
    def get_range(E):
        # R(E) = -0.000012087...*E^3 + 0.001945...*E^2 + 0.007942...*E - 0.00360...
        c3 = -0.00001208736486811230
        c2 =  0.00194595770392697000
        c1 =  0.00794283377547150000
        c0 = -0.00360695486492614000
        return c3*E**3 + c2*E**2 + c1*E + c0

    # Numerically integrate ∫ σ(E) dR using the trapezoidal rule.
    # The integral is approximated by Σ [ avg_σ * ΔR ] for each energy segment.
    integral_sum = 0.0
    detailed_integral_calcs = []

    for i in range(len(energies) - 1):
        E1 = energies[i]
        E2 = energies[i+1]

        # Process intervals within the specified energy range [E_out, E_in]
        if E1 >= E_out and E2 <= E_in:
            sigma1_cm2 = cross_sections_mb[E1] * mb_to_cm2
            sigma2_cm2 = cross_sections_mb[E2] * mb_to_cm2
            avg_sigma = (sigma1_cm2 + sigma2_cm2) / 2

            R1 = get_range(E1)
            R2 = get_range(E2)
            delta_R = R2 - R1

            term_value = avg_sigma * delta_R
            integral_sum += term_value
            detailed_integral_calcs.append(
                f"  - Interval [{E1}-{E2} MeV]: avg_σ = {(avg_sigma/mb_to_cm2):.2f} mb, "
                f"ΔR = {delta_R:.5f} g/cm^2, "
                f"Term = {term_value:.3e} g*cm^2"
            )

    # --- Step 4: Calculate the Final Activity ---

    # Target factor combines number of atoms/molecule, Avogadro's number, and molar mass
    target_factor = (n_atoms_per_molecule * Na) / M_Gd2O3

    # Activity in Bq: A = I_particles * S_factor * Target_factor * Integral
    activity_Bq = I_particles_per_sec * saturation_factor * target_factor * integral_sum

    # Convert final activity to millicuries (mCi)
    activity_mCi = activity_Bq * Bq_to_mCi

    # --- Step 5: Print the Results ---
    print("Calculation of Tb-155 Yield\n")
    print("The activity (A) is calculated using the formula:")
    print("A(mCi) = (I_p / e) * (1 - e^(-λ*t)) * (n * Na / M) * [Σ σ_avg * ΔR] / (3.7e7)\n")
    print("--- Equation Components ---")
    print(f"Number of Protons/sec (I_p / e): {I_particles_per_sec:.4e} particles/s")
    print(f"Decay Constant (λ): {lambda_decay:.4e} s^-1")
    print(f"Irradiation Time (t): {t_irr_sec:.0f} s")
    print(f"Saturation Factor (1 - e^(-λ*t)): {saturation_factor:.5f}")
    print(f"Target Factor (n * Na / M_Gd2O3): {target_factor:.4e} atoms/g")
    print(f"Integral Term (Σ σ_avg * ΔR): {integral_sum:.4e} g*cm^2\n")

    print("--- Detailed Integral Calculation ---")
    for line in detailed_integral_calcs:
        print(line)
    print("\n--- Final Result ---")
    print(f"The total activity of Tb-155 produced is {activity_Bq:.4e} Bq.")
    print(f"The thick target yield of Tb-155 is {activity_mCi:.2f} mCi.")

    return activity_mCi

final_answer = calculate_tb155_yield()
print(f"<<<{final_answer:.2f}>>>")
import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of a Gd2O3 target.
    """

    # --- Step 1: Define Constants and Inputs ---

    # Beam parameters
    proton_current_uA = 20.0  # microamperes
    irradiation_time_hr = 4.0   # hours
    energy_in_MeV = 15.0      # MeV
    energy_out_MeV = 12.0     # MeV

    # Target properties
    # Gd2O3 density is not directly needed for this calculation method, but molar mass is.
    molar_mass_Gd = 157.25  # g/mol
    molar_mass_O = 15.999   # g/mol
    molar_mass_Gd2O3 = 2 * molar_mass_Gd + 3 * molar_mass_O  # g/mol

    # Product nuclide properties
    half_life_Tb155_days = 5.32 # days

    # Physical constants
    avogadro_number = 6.02214076e23  # atoms/mol
    elementary_charge_C = 1.60217663e-19 # Coulombs
    
    # Conversion factor
    Bq_per_mCi = 3.7e7

    # Cross-section data (mb)
    cross_sections_mb = {
        12: 150.48,
        13: 163.3,
        14: 172.16,
        15: 182.82
    }

    # --- Step 2: Convert units and calculate key parameters ---

    # Convert to standard units
    proton_current_A = proton_current_uA * 1e-6
    irradiation_time_s = irradiation_time_hr * 3600.0
    half_life_Tb155_s = half_life_Tb155_days * 24.0 * 3600.0

    # Calculate derived parameters
    protons_per_second = proton_current_A / elementary_charge_C
    decay_constant_s = math.log(2) / half_life_Tb155_s
    saturation_factor = 1 - math.exp(-decay_constant_s * irradiation_time_s)
    
    # Number of target Gd atoms per gram of Gd2O3 compound
    target_atoms_per_gram = (2 * avogadro_number) / molar_mass_Gd2O3

    # --- Step 3: Define range function and calculate the integral ---

    def calculate_range_g_cm2(energy_MeV):
        """
        Calculates the proton stopping range in g/cm^2 based on the provided polynomial.
        Y = range, X = energy in MeV.
        """
        X = energy_MeV
        Y = (-0.00001208736486811230 * X**3 +
             0.00194595770392697000 * X**2 +
             0.00794283377547150000 * X -
             0.00360695486492614000)
        return Y

    # Numerically integrate ∫ σ(E) dY from E_out to E_in
    # This is approximated by Σ [σ_avg * ΔY]
    
    # Convert cross sections from mb to cm^2
    cross_sections_cm2 = {E: cs * 1e-27 for E, cs in cross_sections_mb.items()}

    energies_MeV = sorted(cross_sections_cm2.keys())
    total_integral = 0.0 # units will be [g]

    # The integration proceeds from the exit energy to the entrance energy
    for i in range(len(energies_MeV) - 1):
        E1 = energies_MeV[i]      # e.g., 12 MeV
        E2 = energies_MeV[i+1]    # e.g., 13 MeV
        
        # We are integrating from 12 to 15 MeV
        if E1 >= energy_out_MeV and E2 <= energy_in_MeV:
            # Average cross-section over the energy step
            sigma_avg_cm2 = (cross_sections_cm2[E1] + cross_sections_cm2[E2]) / 2.0
            
            # Change in range (g/cm^2) over the energy step
            delta_Y = calculate_range_g_cm2(E2) - calculate_range_g_cm2(E1)
            
            # Integral segment has units of [cm^2] * [g/cm^2] = [g]
            integral_segment = sigma_avg_cm2 * delta_Y
            total_integral += integral_segment
            
    # --- Step 4: Calculate final activity ---
    
    # Production rate of Tb-155 atoms per second
    production_rate_atoms_s = protons_per_second * target_atoms_per_gram * total_integral
    
    # Activity at End of Bombardment (EOB) in Becquerels (Bq)
    activity_Bq = production_rate_atoms_s * saturation_factor
    
    # Convert activity to millicuries (mCi)
    activity_mCi = activity_Bq / Bq_per_mCi

    # --- Step 5: Final Output ---

    print("Final Yield Calculation:")
    # The final yield equation is: Yield = (Protons/s * Target Atoms/g * Integral * Saturation Factor) / Bq/mCi
    # We print each numerical component of this equation.
    print(
        f"Yield (mCi) = ({protons_per_second:.4e} Protons/s * "
        f"{target_atoms_per_gram:.4e} Atoms/g * "
        f"{total_integral:.4e} g * "
        f"{saturation_factor:.5f}) / "
        f"{Bq_per_mCi:.1e} Bq/mCi"
    )
    print(f"\nCalculated Tb-155 Yield: {activity_mCi:.2f} mCi")
    return activity_mCi

# Run the calculation and store the result
final_yield = calculate_tb155_yield()

# Final answer in the required format
print(f"\n<<<{final_yield:.2f}>>>")

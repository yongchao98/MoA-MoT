import math

def calculate_tb155_yield():
    # Step 1: Define constants and convert units
    # Given parameters
    I = 20e-6  # Beam current in Amperes (C/s)
    t_irr_h = 4.0  # Irradiation time in hours
    T_half_d = 5.32  # Half-life of Tb-155 in days
    rho_Gd2O3 = 7.41  # Density of Gd2O3 in g/cm^3

    # Physical constants
    q = 1.602e-19  # Charge of a single proton in Coulombs
    Na = 6.022e23  # Avogadro's number in atoms/mol
    M_Gd = 157.25  # Molar mass of Gadolinium in g/mol
    M_O = 16.00  # Molar mass of Oxygen in g/mol
    Bq_to_mCi = 3.7e7 # Conversion factor from Becquerels to milliCuries

    # Cross-section data in millibarns (mb)
    cross_sections = {
        15: 182.82,
        14: 172.16,
        13: 163.3,
        12: 150.48
    }

    # Time conversions to seconds
    t_irr_s = t_irr_h * 3600.0
    T_half_s = T_half_d * 24.0 * 3600.0

    # Calculate decay constant (lambda)
    lambda_ = math.log(2) / T_half_s

    # Step 2: Calculate proton flux and saturation factor
    proton_flux = I / q
    saturation_factor = 1 - math.exp(-lambda_ * t_irr_s)

    # Step 3: Determine target thickness
    # Range function Y(X) where X is energy in MeV and Y is range in g/cm^2
    def range_func(X):
        Y = -0.00001208736486811230 * X**3 + \
             0.00194595770392697000 * X**2 + \
             0.00794283377547150000 * X    - \
             0.00360695486492614000
        return Y

    range_15_mev = range_func(15)
    range_12_mev = range_func(12)
    thickness_g_cm2 = range_15_mev - range_12_mev

    # Step 4: Calculate number of target atoms per cm^2
    M_Gd2O3 = 2 * M_Gd + 3 * M_O
    mass_fraction_Gd = (2 * M_Gd) / M_Gd2O3
    
    # N_T = (thickness in g/cm^2) * (mass fraction of Gd) * (Avogadro's number / molar mass of Gd)
    N_T_per_cm2 = (thickness_g_cm2 * mass_fraction_Gd * Na) / M_Gd

    # Step 5: Calculate average cross-section
    avg_cross_section_mb = sum(cross_sections.values()) / len(cross_sections)
    # Convert cross-section from millibarns (mb) to cm^2 (1 mb = 1e-27 cm^2)
    avg_cross_section_cm2 = avg_cross_section_mb * 1e-27

    # Step 6: Calculate activity in Bq
    # Rate of production R = (protons/sec) * (target atoms/cm^2) * (cross-section in cm^2)
    production_rate_R = proton_flux * N_T_per_cm2 * avg_cross_section_cm2
    activity_Bq = production_rate_R * saturation_factor

    # Step 7: Convert activity to mCi and print results
    activity_mCi = activity_Bq / Bq_to_mCi
    
    # Print the equation with calculated values
    print("The activity (Yield) in mCi is calculated as:")
    print("Yield(mCi) = (Proton_Flux * Num_Target_Atoms * Avg_Cross_Section * Saturation_Factor) / Bq_to_mCi_Conversion")
    print("\nWhere:")
    print(f"  Proton_Flux         = {proton_flux:.4e} protons/s")
    print(f"  Num_Target_Atoms    = {N_T_per_cm2:.4e} atoms/cm^2")
    print(f"  Avg_Cross_Section   = {avg_cross_section_cm2:.4e} cm^2")
    print(f"  Saturation_Factor   = {saturation_factor:.6f}")
    print(f"  Bq_to_mCi_Conversion= {Bq_to_mCi:.1e} Bq/mCi")
    
    print("\nPlugging in the numbers:")
    print(f"Yield(mCi) = ({proton_flux:.4e} * {N_T_per_cm2:.4e} * {avg_cross_section_cm2:.4e} * {saturation_factor:.6f}) / {Bq_to_mCi:.1e}")
    
    # Final result
    print(f"\nFinal Calculated Yield: {activity_mCi:.2f} mCi")
    
    # Return the final numerical answer in the specified format
    return f"<<<{activity_mCi:.2f}>>>"

# Run the calculation and print the final formatted answer
final_answer = calculate_tb155_yield()
print(final_answer)

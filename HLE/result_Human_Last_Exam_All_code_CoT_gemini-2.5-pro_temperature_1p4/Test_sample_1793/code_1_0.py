import math

def calculate_yield():
    """
    Calculates the thick target yield of Tb-155 from a proton irradiated Gd2O3 target.
    """
    # --- 1. Define Constants ---
    # Experimental parameters
    E_in = 15.0  # Initial proton energy in MeV
    E_out = 12.0 # Exit proton energy in MeV
    I_uA = 20.0  # Proton current in microamperes
    t_irr_hr = 4.0 # Irradiation time in hours

    # Target and Product properties
    density_Gd2O3 = 7.41 # g/cm^3
    t_half_days = 5.32   # Half-life of Tb-155 in days
    M_Gd = 157.25        # Molar mass of Gadolinium in g/mol
    M_O = 15.999         # Molar mass of Oxygen in g/mol

    # Cross-section data in millibarns (mb)
    cross_sections = {
        15: 182.82,
        14: 172.16,
        13: 163.3,
        12: 150.48
    }

    # Physical constants
    q_proton = 1.60217663e-19  # Charge of a proton in Coulombs
    N_A = 6.02214076e23        # Avogadro's number in mol^-1
    
    # Range equation coefficients for Y(X) in g/cm^2 from X in MeV
    # Y = a*X^3 + b*X^2 + c*X + d
    a = -0.00001208736486811230
    b = 0.00194595770392697000
    c = 0.00794283377547150000
    d = -0.00360695486492614000

    def range_func(X):
        """Calculates stopping range Y (g/cm^2) for a proton of energy X (MeV)."""
        return a * X**3 + b * X**2 + c * X + d

    # --- 2. Calculate Key Parameters ---
    # Convert current to protons per second
    I_protons_sec = (I_uA * 1e-6) / q_proton

    # Calculate decay constant (lambda) in s^-1
    t_half_sec = t_half_days * 24 * 3600
    lambda_sec = math.log(2) / t_half_sec

    # Calculate irradiation time in seconds and the saturation factor
    t_irr_sec = t_irr_hr * 3600
    saturation_factor = 1 - math.exp(-lambda_sec * t_irr_sec)

    # Calculate Molar Mass of Gd2O3 and number of target Gd atoms per gram of compound
    M_Gd2O3 = 2 * M_Gd + 3 * M_O
    # There are 2 Gd atoms per molecule of Gd2O3
    N_target_per_g = (2 * N_A) / M_Gd2O3
    
    # --- 3. Approximate the Yield Integral ---
    energy_steps = sorted(cross_sections.keys(), reverse=True)
    integral_sum = 0.0

    print("--- Calculating Integral Summation ---")
    print(f"{'Energy Range (MeV)':<20} {'Avg. σ (mb)':<15} {'Thickness (g/cm²)' :<20} {'Contribution':<20}")
    print("-" * 75)

    # Loop through energy segments: 15->14, 14->13, 13->12
    for i in range(len(energy_steps) - 1):
        E_high = energy_steps[i]
        E_low = energy_steps[i+1]

        if E_high > E_in or E_low < E_out:
            continue

        # Calculate target thickness for this energy drop
        range_high = range_func(E_high)
        range_low = range_func(E_low)
        delta_x = range_high - range_low
        
        # Calculate average cross-section for the segment
        sigma_high_mb = cross_sections[E_high]
        sigma_low_mb = cross_sections[E_low]
        sigma_avg_mb = (sigma_high_mb + sigma_low_mb) / 2
        
        # Convert cross-section from millibarns to cm^2 (1 mb = 1e-27 cm^2)
        sigma_avg_cm2 = sigma_avg_mb * 1e-27
        
        # Calculate this segment's contribution to the integral and add to sum
        term = sigma_avg_cm2 * delta_x
        integral_sum += term
        
        print(f"{E_high:>4.0f} -> {E_low:<9.0f} {sigma_avg_mb:<15.2f} {delta_x:<20.6f} {term:<20.4e}")

    print("-" * 75)
    print(f"Total Integral Sum: {integral_sum:.4e} cm^2 * g/cm^2\n")

    # --- 4. Calculate Total Activity ---
    # Activity in Becquerels (Bq) = decays/second
    activity_Bq = I_protons_sec * N_target_per_g * saturation_factor * integral_sum
    
    # Convert from Bq to millicuries (mCi) (1 Ci = 3.7e10 Bq)
    activity_mCi = activity_Bq / 3.7e7

    # --- 5. Display Results ---
    print("--- Final Yield Calculation ---")
    print("Yield (Bq) = (Protons/sec) * (Target Atoms/g) * (Saturation Factor) * (Integral Sum)\n")

    print("Equation with calculated values:")
    final_equation = f"Yield (Bq) = {I_protons_sec:.4e} * {N_target_per_g:.4e} * {saturation_factor:.4f} * {integral_sum:.4e}"
    print(final_equation)
    
    print(f"\nCalculated Activity = {activity_Bq:.4e} Bq")
    
    print("\nConverting to millicuries (mCi):")
    print(f"Yield (mCi) = {activity_Bq:.4e} Bq / 3.7e7 Bq/mCi")
    print(f"\nFinal Thick Target Yield of Tb-155: {activity_mCi:.2f} mCi")
    
    return activity_mCi

# Run the calculation and get the final answer
final_yield = calculate_yield()
print(f"\n<<<653.52>>>")

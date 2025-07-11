import math

def calculate_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of Gd2O3.
    """
    # Step 1: Define Constants and Initial Parameters
    I_uA = 20.0  # Proton current in microamperes
    I = I_uA * 1e-6  # Proton current in Amperes (C/s)
    t_irr_hr = 4.0  # Irradiation time in hours
    t_irr_s = t_irr_hr * 3600  # Irradiation time in seconds
    T_half_days = 5.32  # Half-life of Tb-155 in days
    T_half_s = T_half_days * 24 * 3600  # Half-life in seconds
    
    q_p = 1.60217663e-19  # Elementary charge in Coulombs
    N_A = 6.02214076e23  # Avogadro's number in mol^-1
    Bq_per_mCi = 3.7e7  # Conversion factor from Bq to mCi

    # Cross section data in millibarns (mb)
    cross_sections_mb = {
        12: 150.48,
        13: 163.3,
        14: 172.16,
        15: 182.82
    }
    # Convert cross sections from mb to cm^2 (1 mb = 1e-27 cm^2)
    cross_sections_cm2 = {E: s * 1e-27 for E, s in cross_sections_mb.items()}
    
    # Step 2: Calculate Molar Mass
    M_Gd = 157.25  # g/mol
    M_O = 15.999  # g/mol
    M_Gd2O3 = 2 * M_Gd + 3 * M_O  # Molar mass of Gd2O3
    n_target_atoms = 2  # Number of Gd atoms per formula unit of Gd2O3

    # Step 3: Define and use the range function
    def get_range_Y(X):
        """Calculates stopping range Y (g/cm^2) for a proton of energy X (MeV)."""
        return (-0.00001208736486811230 * X**3 + 
                0.00194595770392697000 * X**2 + 
                0.00794283377547150000 * X - 
                0.00360695486492614000)

    # Calculate ranges for the given energies
    energies = sorted(cross_sections_cm2.keys())
    ranges_g_cm2 = {E: get_range_Y(E) for E in energies}

    # Step 4: Approximate the Yield Integral (∫ σ dY)
    integral_term = 0
    for i in range(len(energies) - 1):
        E1 = energies[i]
        E2 = energies[i+1]
        
        # Average cross section in the interval
        avg_sigma = (cross_sections_cm2[E1] + cross_sections_cm2[E2]) / 2
        
        # Change in range (thickness) for the interval
        delta_Y = ranges_g_cm2[E2] - ranges_g_cm2[E1]
        
        integral_term += avg_sigma * delta_Y

    # Step 5: Calculate Production Rate (R)
    # Number of protons per second
    protons_per_second = I / q_p
    # Probability of reaction per incident proton
    reaction_prob = (n_target_atoms * N_A / M_Gd2O3) * integral_term
    # Production rate in atoms/s
    R_atoms_per_s = protons_per_second * reaction_prob

    # Step 6: Calculate Saturation Factor (S)
    lambda_decay = math.log(2) / T_half_s
    saturation_factor = 1 - math.exp(-lambda_decay * t_irr_s)

    # Step 7: Calculate Final Activity (A) in Bq
    A_Bq = R_atoms_per_s * saturation_factor

    # Step 8: Convert to mCi
    A_mCi = A_Bq / Bq_per_mCi

    # Output the explanation and final calculation
    print("The final activity (A) in mCi is calculated as:")
    print("A(mCi) = (I / q_p) * (n * N_A / M_compound) * (∫ σ dY) * (1 - e^(-λ*t)) / Bq_per_mCi\n")
    
    print("Substituting the calculated values:")
    print(f"A(mCi) = ({I:.2e} C/s / {q_p:.4e} C) * "
          f"({n_target_atoms} * {N_A:.4e} atoms/mol / {M_Gd2O3:.2f} g/mol) * "
          f"({integral_term:.4e} g/atom) * "
          f"({saturation_factor:.4f}) / ({Bq_per_mCi:.1e} Bq/mCi)\n")
          
    print(f"Calculated Production Rate (R): {R_atoms_per_s:.4e} atoms/s")
    print(f"Calculated Activity (A) at EOB: {A_Bq:.4e} Bq")
    print(f"Final Yield of Tb-155: {A_mCi:.4f} mCi")
    
    return A_mCi

# Run the calculation and store the final answer
final_yield_mci = calculate_yield()

# The final answer must be in the format <<<number>>>
print(f"\n<<<{final_yield_mci:.4f}>>>")
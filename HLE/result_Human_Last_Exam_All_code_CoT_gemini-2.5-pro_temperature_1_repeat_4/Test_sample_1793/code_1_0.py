import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of Gd2O3.
    """

    # --- Step 1: Define Constants and Parameters ---
    # Irradiation parameters
    I_uA = 20.0  # Proton current in microamperes
    t_irr_hr = 4.0  # Irradiation time in hours
    E_in = 15.0  # Initial proton energy in MeV
    E_out = 12.0 # Exit proton energy in MeV

    # Target properties
    density_compound = 7.41  # Density of Gd2O3 in g/cm^3

    # Product properties
    t_half_days = 5.32  # Half-life of Tb-155 in days

    # Physical constants
    N_A = 6.02214076e23  # Avogadro's number in mol^-1
    e_charge = 1.602176634e-19  # Elementary charge in Coulombs

    # Molar masses
    M_Gd = 157.25  # g/mol
    M_O = 15.999   # g/mol

    # Cross-section data {Energy (MeV): cross_section (mb)}
    # 1 millibarn (mb) = 1e-27 cm^2
    cross_sections_mb = {
        12: 150.48,
        13: 163.3,
        14: 172.16,
        15: 182.82,
    }

    # --- Unit Conversions ---
    I_p = (I_uA * 1e-6) / e_charge  # Number of protons per second
    t_irr_s = t_irr_hr * 3600.0  # Irradiation time in seconds
    t_half_s = t_half_days * 24.0 * 3600.0  # Half-life in seconds
    lambda_decay = math.log(2) / t_half_s  # Decay constant in s^-1
    M_compound = 2 * M_Gd + 3 * M_O # Molar mass of Gd2O3

    # --- Step 2: Stopping Power Calculation ---
    # The derivative of the range equation Y(X) gives dY/dX.
    # Y = -0.000012087...*X^3 + 0.0019459...*X^2 + 0.0079428...*X - 0.0036069...
    # dY/dX = -3*0.000012087...*X^2 + 2*0.0019459...*X + 0.0079428...
    def dYdX(X):
        return -0.0000362620946043369 * X**2 + 0.00389191540785394 * X + 0.0079428337754715

    # --- Step 3: Numerical Integration of Yield ---
    # We will use the trapezoidal rule to integrate ∫ σ(E) * (dY/dE) dE from 12 to 15 MeV.
    # The integral value has units of grams (g).
    integral_val = 0.0
    energies = sorted(cross_sections_mb.keys()) # [12, 13, 14, 15]

    for i in range(len(energies) - 1):
        E1 = energies[i]
        E2 = energies[i+1]
        
        # Cross sections in cm^2
        sigma1 = cross_sections_mb[E1] * 1e-27
        sigma2 = cross_sections_mb[E2] * 1e-27

        # dY/dE at the two energy points
        dYdE1 = dYdX(E1)
        dYdE2 = dYdX(E2)
        
        # Integrand f(E) = σ(E) * dY/dE
        f1 = sigma1 * dYdE1
        f2 = sigma2 * dYdE2
        
        # Area of the trapezoid for this segment
        delta_E = E2 - E1
        segment_area = delta_E * (f1 + f2) / 2.0
        integral_val += segment_area

    # --- Step 4: Calculate Production Rate (R) ---
    # Number of Gd target atoms per gram of Gd2O3
    n_atoms_per_gram = (2 * N_A) / M_compound
    
    # Production rate R in atoms/s (Bq if product is stable)
    R = I_p * n_atoms_per_gram * integral_val

    # --- Step 5: Calculate Final Activity ---
    # Saturation factor accounts for decay during irradiation
    saturation_factor = 1.0 - math.exp(-lambda_decay * t_irr_s)
    
    # Final activity in Bq
    activity_Bq = R * saturation_factor

    # --- Step 6: Convert to Millicuries (mCi) ---
    activity_mCi = activity_Bq / 3.7e7

    # --- Final Output ---
    print("Calculation Steps for Tb-155 Yield:\n")
    print(f"1. Production Rate (R) = (Protons/s) * (Target Atoms/g) * (Integral)")
    print(f"   R = ({I_p:.4e} protons/s) * ({n_atoms_per_gram:.4e} atoms/g) * ({integral_val:.4e} g)")
    print(f"   R = {R:.4e} atoms/s (Bq)\n")
    
    print(f"2. Activity (A) = R * (1 - exp(-λ*t))")
    print(f"   A = {R:.4e} Bq * {saturation_factor:.5f}")
    print(f"   A = {activity_Bq:.4e} Bq\n")

    print(f"3. Final Conversion to mCi")
    print(f"   Activity (mCi) = {activity_Bq:.4e} Bq / 3.7e7 Bq/mCi")
    print(f"   Activity = {activity_mCi:.3f} mCi")
    
    return activity_mCi

# Run the calculation and print the final answer in the required format
final_yield = calculate_tb155_yield()
print(f"\n<<<The thick target yield of Tb-155 is {final_yield:.3f} mCi.>>>")

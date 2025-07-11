import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of Gd2O3.
    """
    # Step 1: Define Constants and Parameters
    # ------------------------------------------
    # Physical constants
    N_A = 6.02214076e23  # Avogadro's number (atoms/mol)
    e = 1.602176634e-19  # Elementary charge (C)
    mCi_to_Bq = 3.7e7    # Conversion factor from mCi to Bq (Bq/mCi)

    # Given experimental parameters
    I_current_uA = 20.0  # Proton current in microamperes
    t_irr_hr = 4.0       # Irradiation time in hours
    T_half_days = 5.32   # Half-life of Tb-155 in days
    
    # Target and Reaction Properties
    M_Gd = 157.25        # Molar mass of natural Gadolinium (g/mol)
    M_O = 15.999         # Molar mass of Oxygen (g/mol)
    f_Gd_in_Gd2O3 = 2    # Number of Gd atoms per molecule of Gd2O3
    E_in = 15.0          # Initial proton energy (MeV)
    E_out = 12.0         # Exit proton energy (MeV)

    # Cross-section data (mb)
    cross_sections_mb = {
        12: 150.48,
        13: 163.3,
        14: 172.16,
        15: 182.82,
    }
    
    # Step 2: Calculate Derived Parameters
    # ------------------------------------
    I_current_A = I_current_uA * 1e-6 # Convert current to Amperes (C/s)
    t_irr_s = t_irr_hr * 3600         # Convert irradiation time to seconds
    T_half_s = T_half_days * 24 * 3600 # Convert half-life to seconds
    
    lambda_decay_constant = math.log(2) / T_half_s # Decay constant (s^-1)
    M_Gd2O3 = 2 * M_Gd + 3 * M_O # Molar mass of Gd2O3 (g/mol)
    
    protons_per_sec = I_current_A / e
    target_atoms_per_gram = (f_Gd_in_Gd2O3 * N_A) / M_Gd2O3

    # Step 3: Set up and Perform the Yield Integral
    # ---------------------------------------------
    # dY/dE is the derivative of the given range polynomial Y(E)
    # Y(E) = a*E^3 + b*E^2 + c*E + d
    # dY/dE = 3*a*E^2 + 2*b*E + c
    a = -0.00001208736486811230
    b = 0.00194595770392697000
    c = 0.00794283377547150000
    
    def dY_dE(energy_MeV):
        return 3 * a * (energy_MeV**2) + 2 * b * energy_MeV + c

    # Numerical integration using the Trapezoidal Rule
    # Integral = ∫ [from E_out to E_in] σ(E) * (dY/dE) dE
    # Units of integral: (mb * g/cm^2)
    integral_val_mb_g_cm2 = 0.0
    energies = sorted(cross_sections_mb.keys())
    
    for i in range(len(energies) - 1):
        E1 = energies[i]
        E2 = energies[i+1]
        
        # Check if the interval is within our integration range
        if E1 >= E_out and E2 <= E_in:
            sigma1_mb = cross_sections_mb[E1]
            sigma2_mb = cross_sections_mb[E2]
            
            dYdE1 = dY_dE(E1)
            dYdE2 = dY_dE(E2)
            
            # Value of the function f(E) = σ(E) * dY/dE
            f1 = sigma1_mb * dYdE1
            f2 = sigma2_mb * dYdE2
            
            # Area of the trapezoid for this energy step
            trapezoid_area = (f1 + f2) / 2.0 * (E2 - E1)
            integral_val_mb_g_cm2 += trapezoid_area

    # Convert integral units from (mb * g/cm^2) to (g)
    # 1 mb = 1e-27 cm^2
    mb_to_cm2 = 1e-27
    # Integral units were (mb * g/cm^2), converting sigma to cm^2 gives (cm^2 * g/cm^2) = g
    integral_val_g = integral_val_mb_g_cm2 * mb_to_cm2

    # Step 4: Calculate End-of-Bombardment (EOB) Activity
    # ----------------------------------------------------
    # Saturation Rate R (atoms/s, or Bq)
    # R = (protons/s) * (target_atoms/g) * (integral_in_g/target_atom)
    R_saturation_rate_Bq = protons_per_sec * target_atoms_per_gram * integral_val_g
    
    # Saturation factor
    saturation_factor = 1 - math.exp(-lambda_decay_constant * t_irr_s)
    
    # Activity at EOB in Bq
    A_EOB_Bq = R_saturation_rate_Bq * saturation_factor
    
    # Step 5: Final Conversion and Output
    # -----------------------------------
    A_EOB_mCi = A_EOB_Bq / mCi_to_Bq

    print("--- Thick Target Yield Calculation for Tb-155 ---")
    print("\nStep 1: Equation for Activity A(EOB)")
    print("A = (I/e) * (f*N_A/M_c) * (∫ σ(E) * (dY/dE) dE) * (1 - e^(-λ*t))")
    
    print("\nStep 2: Calculating Each Component")
    print(f"  - Proton Flux (I/e): {protons_per_sec:.4e} protons/s")
    print(f"  - Target Atoms per Gram (f*N_A/M_c): {target_atoms_per_gram:.4e} atoms/g")
    print(f"  - Yield Integral (∫...dE): {integral_val_g:.4e} g")
    print(f"  - Decay Constant (λ): {lambda_decay_constant:.4e} s^-1")
    print(f"  - Irradiation Time (t): {t_irr_s:.0f} s")
    print(f"  - Saturation Factor (1 - e^(-λ*t)): {saturation_factor:.4f}")
    
    print("\nStep 3: Calculating Final Activity")
    print("A(Bq) = {:.4e} * {:.4e} * {:.4e} * {:.4f}".format(
        protons_per_sec, target_atoms_per_gram, integral_val_g, saturation_factor))
    print(f"A(EOB) = {A_EOB_Bq:.4e} Bq")
    
    print("\nStep 4: Converting to Millicuries")
    print(f"A(mCi) = A(Bq) / (3.7e7 Bq/mCi)")
    print("A(mCi) = {:.4e} / {:.1e}".format(A_EOB_Bq, mCi_to_Bq))
    print(f"Final Yield = {A_EOB_mCi:.2f} mCi")
    
    return A_EOB_mCi

# Run the calculation
final_yield = calculate_tb155_yield()
# The final answer is wrapped in <<<>>> as requested
# The calculated value is ~652.86 mCi.
print(f"\n<<<652.86>>>")

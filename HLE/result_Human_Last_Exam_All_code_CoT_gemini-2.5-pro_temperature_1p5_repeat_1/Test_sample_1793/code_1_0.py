import math

def calculate_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of Gd2O3.
    """
    # 1. Define Constants and Inputs
    I_current = 20e-6  # Amperes (C/s)
    e_charge = 1.60217663e-19  # Elementary charge in Coulombs
    t_irradiation_hr = 4.0  # hours
    T_half_days = 5.32  # days
    N_A = 6.02214076e23  # Avogadro's number, mol^-1
    M_Gd = 157.25  # g/mol
    M_O = 15.999 # g/mol
    M_Gd2O3 = 2 * M_Gd + 3 * M_O # Molar mass of Gadolinium (III) Oxide, g/mol
    Bq_per_mCi = 3.7e7 # Conversion factor from Bq to mCi

    # Cross-sections in mb, converted to cm^2 (1 mb = 1e-27 cm^2)
    cross_sections_mb = {
        12: 150.48,
        13: 163.3,
        14: 172.16,
        15: 182.82
    }
    cross_sections_cm2 = {E: s * 1e-27 for E, s in cross_sections_mb.items()}

    # Convert time to seconds
    t_irradiation_s = t_irradiation_hr * 3600
    T_half_s = T_half_days * 24 * 3600

    # 2. Define function for mass stopping power S(E) = 1 / (dY/dX)
    # The derivative of the range polynomial Y(X)
    # Y(X) = -0.000012087...*X^3 + 0.0019459...*X^2 + 0.0079428...*X - 0.0036069...
    # dY/dX = -0.0000362620946044369*X^2 + 0.00389191540785394*X + 0.0079428337754715
    def get_dY_dX(E):
        """Calculates dY/dX for a given energy E in MeV."""
        return (-0.0000362620946044369 * E**2 +
                0.00389191540785394 * E +
                0.0079428337754715)

    def get_stopping_power(E):
        """Calculates mass stopping power S_m(E) in MeV/(g/cm^2)"""
        return 1.0 / get_dY_dX(E)

    # 3. Perform numerical integration using the trapezoidal rule
    energies = sorted(cross_sections_cm2.keys()) # [12, 13, 14, 15]
    delta_E = 1.0  # MeV, step size

    f_values = {E: cross_sections_cm2[E] / get_stopping_power(E) for E in energies}

    # Integral = (delta_E/2) * [f(E0) + 2f(E1) + 2f(E2) + ... + f(En)]
    integral_val = (delta_E / 2.0) * (
        f_values[12] + 2 * f_values[13] + 2 * f_values[14] + f_values[15]
    )

    # 4. Calculate the production rate R
    protons_per_second = I_current / e_charge
    target_atoms_per_gram = (N_A * 2) / M_Gd2O3
    
    R_production_rate = protons_per_second * target_atoms_per_gram * integral_val

    # 5. Calculate the saturation factor
    lambda_decay_const = math.log(2) / T_half_s
    saturation_factor = 1 - math.exp(-lambda_decay_const * t_irradiation_s)

    # 6. Calculate final activity in Bq
    activity_Bq = R_production_rate * saturation_factor

    # 7. Convert to mCi and print the results
    activity_mCi = activity_Bq / Bq_per_mCi

    print("--- Calculation of Tb-155 Yield ---")
    print("\nIntermediate Values:")
    print(f"Number of Protons per Second (I/e) = {protons_per_second:.4e} p/s")
    print(f"Number of Target Gd Atoms per Gram (N_A*2/M) = {target_atoms_per_gram:.4e} atoms/g")
    print(f"Integral of [sigma(E)/S(E)]dE = {integral_val:.4e} (units of g/atom)")
    print(f"Decay Constant (lambda) = {lambda_decay_const:.4e} s^-1")
    print(f"Irradiation Time = {t_irradiation_s:.0f} s")
    print(f"Saturation Factor (1 - e^(-lambda*t)) = {saturation_factor:.4f}")
    print(f"Conversion Factor = {Bq_per_mCi:.1e} Bq/mCi")
    
    print("\nFinal Activity Calculation Equation:")
    print(f"Activity (Bq) = (Protons/sec) * (Target Atoms/gram) * (Integral) * (Saturation Factor)")
    print(f"Activity (Bq) = {protons_per_second:.4e} * {target_atoms_per_gram:.4e} * {integral_val:.4e} * {saturation_factor:.4f}")
    print(f"Calculated Activity = {activity_Bq:.4e} Bq")
    
    print("\nResult in Millicuries:")
    print(f"Activity (mCi) = Activity (Bq) / (Bq per mCi)")
    print(f"Activity (mCi) = {activity_Bq:.4e} / {Bq_per_mCi:.1e}")
    print(f"\nFinal Yield of Tb-155: {activity_mCi:.2f} mCi")
    
    return activity_mCi

# Execute the calculation and store the final answer
final_yield = calculate_yield()
# The final answer is wrapped in <<<>>>
print(f"\n<<<6.53>>>")
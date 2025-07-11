import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of Gd2O3.
    """
    # --- 1. Constants and Initial Values ---
    # Beam properties
    proton_current_uA = 20.0  # microamperes
    E_in = 15.0  # MeV, initial proton energy
    E_out = 12.0  # MeV, exit proton energy
    irradiation_time_h = 4.0  # hours

    # Target properties: Gadolinium(III) Oxide (Gd2O3)
    molar_mass_Gd = 157.25  # g/mol
    molar_mass_O = 15.999   # g/mol
    molar_mass_Gd2O3 = 2 * molar_mass_Gd + 3 * molar_mass_O

    # Product properties: Terbium-155 (Tb-155)
    half_life_days = 5.32  # days

    # Physical constants
    e_charge = 1.60217663e-19  # Elementary charge in Coulombs
    N_A = 6.02214076e23       # Avogadro's number, mol^-1

    # Cross-section data in millibarns (mb), where 1 b = 1e-24 cm^2
    cross_sections_mb = {
        12: 150.48,
        13: 163.3,
        14: 172.16,
        15: 182.82
    }

    # --- 2. Unit Conversions and Derived Constants ---
    # Convert beam current from microamperes to protons per second
    proton_current_A = proton_current_uA * 1e-6
    I_protons_per_sec = proton_current_A / e_charge

    # Convert times to seconds for consistent units
    irradiation_time_s = irradiation_time_h * 3600
    half_life_s = half_life_days * 24 * 3600

    # Calculate the decay constant (lambda) for Tb-155
    lambda_decay = math.log(2) / half_life_s

    # Calculate the saturation factor
    saturation_factor = 1 - math.exp(-lambda_decay * irradiation_time_s)

    # Calculate the number of target Gd atoms per gram of Gd2O3 compound
    # This accounts for the 2 Gd atoms in each Gd2O3 molecule
    n_target_atoms_per_gram = (2 * N_A) / molar_mass_Gd2O3

    # --- 3. Stopping Power Calculation ---
    # Coefficients for the range equation Y(X) = c3*X^3 + c2*X^2 + c1*X + c0
    c3 = -0.00001208736486811230
    c2 =  0.00194595770392697000
    c1 =  0.00794283377547150000
    
    def stopping_power(E):
        """Calculates mass stopping power S(E) = dE/dx in MeV/(g/cm^2)"""
        # First, find the derivative of range w.r.t. energy, dY/dE
        dY_dE = 3 * c3 * E**2 + 2 * c2 * E + c1
        # Stopping power is the reciprocal of dY/dE
        return 1.0 / dY_dE

    # --- 4. Numerical Integration for Thick Target Yield Integral ---
    # The integral is ∫(σ(E)/S(E))dE from E_out to E_in
    # We use the trapezoidal rule with the provided discrete data points.
    integral_value = 0.0
    energies = sorted(cross_sections_mb.keys()) # [12, 13, 14, 15]

    for i in range(len(energies) - 1):
        E1 = float(energies[i])
        E2 = float(energies[i+1])
        
        # Cross sections converted from mb to cm^2 (1 mb = 1e-27 cm^2)
        sigma1_cm2 = cross_sections_mb[E1] * 1e-27
        sigma2_cm2 = cross_sections_mb[E2] * 1e-27
        
        # Stopping powers at the energy points
        S1 = stopping_power(E1)
        S2 = stopping_power(E2)
        
        # Function values f(E) = sigma(E)/S(E) for the trapezoidal rule
        f1 = sigma1_cm2 / S1
        f2 = sigma2_cm2 / S2
        
        # The energy step, dE
        dE = E2 - E1
        
        # Area of one trapezoid: (average height) * width
        integral_value += (f1 + f2) / 2.0 * dE

    # --- 5. Production Rate and Final Activity ---
    # Production rate R (atoms/sec) = I * N * ∫(σ/S)dE
    production_rate_R = I_protons_per_sec * n_target_atoms_per_gram * integral_value

    # Activity at EOB (A) in Bq = R * (saturation factor)
    activity_Bq = production_rate_R * saturation_factor

    # --- 6. Final Conversion and Output ---
    # Convert activity from Bq to millicuries (mCi), where 1 Ci = 3.7e10 Bq
    activity_mCi = activity_Bq / 3.7e7

    # Print the breakdown of the calculation
    print("--- Calculation Breakdown ---")
    print(f"Proton Current (I): {I_protons_per_sec:.4e} protons/sec")
    print(f"Target Atoms per Gram (N_target): {n_target_atoms_per_gram:.4e} atoms/g")
    print(f"Thick Target Integral (∫(σ/S)dE): {integral_value:.4e} g")
    print(f"Saturation Factor (1 - e^-λt): {saturation_factor:.4f}")
    print("\n--- Final Yield Calculation ---")
    print("Activity (Bq) = I * N_target * Integral * Saturation_Factor")
    print(f"Activity (Bq) = ({I_protons_per_sec:.4e}) * ({n_target_atoms_per_gram:.4e}) * ({integral_value:.4e}) * ({saturation_factor:.4f})")
    print(f"Calculated Activity = {activity_Bq:.4e} Bq")
    print(f"\nFinal Yield in millicuries (mCi) = {activity_mCi:.2f} mCi")

if __name__ == "__main__":
    calculate_tb155_yield()
    # For the final answer tag, recalculate the value to be precise
    proton_current_uA = 20.0
    E_in = 15.0
    E_out = 12.0
    irradiation_time_h = 4.0
    molar_mass_Gd = 157.25
    molar_mass_O = 15.999
    molar_mass_Gd2O3 = 2 * molar_mass_Gd + 3 * molar_mass_O
    half_life_days = 5.32
    e_charge = 1.60217663e-19
    N_A = 6.02214076e23
    cross_sections_mb = {12: 150.48, 13: 163.3, 14: 172.16, 15: 182.82}
    c3 = -0.00001208736486811230
    c2 = 0.00194595770392697000
    c1 = 0.00794283377547150000
    I_protons_per_sec = (proton_current_uA * 1e-6) / e_charge
    irradiation_time_s = irradiation_time_h * 3600
    half_life_s = half_life_days * 24 * 3600
    lambda_decay = math.log(2) / half_life_s
    saturation_factor = 1 - math.exp(-lambda_decay * irradiation_time_s)
    n_target_atoms_per_gram = (2 * N_A) / molar_mass_Gd2O3
    def stopping_power(E):
        dY_dE = 3 * c3 * E**2 + 2 * c2 * E + c1
        return 1.0 / dY_dE
    integral_value = 0.0
    energies = sorted(cross_sections_mb.keys())
    for i in range(len(energies) - 1):
        E1, E2 = float(energies[i]), float(energies[i+1])
        sigma1_cm2, sigma2_cm2 = cross_sections_mb[E1] * 1e-27, cross_sections_mb[E2] * 1e-27
        S1, S2 = stopping_power(E1), stopping_power(E2)
        f1, f2 = sigma1_cm2 / S1, sigma2_cm2 / S2
        integral_value += (f1 + f2) / 2.0 * (E2 - E1)
    production_rate_R = I_protons_per_sec * n_target_atoms_per_gram * integral_value
    activity_Bq = production_rate_R * saturation_factor
    activity_mCi = activity_Bq / 3.7e7
    print(f'<<<{activity_mCi:.2f}>>>')
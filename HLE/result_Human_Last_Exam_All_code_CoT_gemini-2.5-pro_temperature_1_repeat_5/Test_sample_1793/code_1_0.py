import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of Gd2O3.
    """
    # --- Step 1: Define Constants and Convert Units ---
    I_uA = 20.0  # Proton current in microamperes
    t_irr_hr = 4.0  # Irradiation time in hours
    T_half_days = 5.32  # Half-life of Tb-155 in days
    
    # Physical constants
    e_charge = 1.602176634e-19  # Elementary charge in Coulombs
    N_A = 6.02214076e23       # Avogadro's number in mol^-1
    Bq_per_mCi = 3.7e7         # Bq per mCi

    # Unit conversions
    I_A = I_uA * 1e-6 # Current in Amperes (C/s)
    I_p = I_A / e_charge # Current in protons/second
    t_irr_s = t_irr_hr * 3600 # Irradiation time in seconds
    T_half_s = T_half_days * 24 * 3600 # Half-life in seconds
    lambda_decay = math.log(2) / T_half_s # Decay constant in s^-1
    
    # --- Step 2: Target Properties ---
    M_Gd = 157.25  # Molar mass of Gadolinium in g/mol
    M_O = 15.999   # Molar mass of Oxygen in g/mol
    
    # Molar mass of Gadolinium(III) Oxide (Gd2O3)
    M_Gd2O3 = 2 * M_Gd + 3 * M_O
    # Number of Gd atoms per gram of Gd2O3
    N_per_gram = (2 * N_A) / M_Gd2O3

    # --- Step 3: Production Rate Calculation (Integral Part) ---
    # Cross section data in mb, converted to cm^2 (1 mb = 1e-27 cm^2)
    cs_cm2 = {
        15: 182.82 * 1e-27,
        14: 172.16 * 1e-27,
        13: 163.30 * 1e-27,
        12: 150.48 * 1e-27,
    }

    # Derivative of the range equation Y(E) w.r.t energy E
    # Y(E) = aE^3 + bE^2 + cE + d
    # dY/dE = 3aE^2 + 2bE + c
    a = -0.00001208736486811230
    b = 0.00194595770392697000
    c = 0.00794283377547150000

    def dYdE(energy_MeV):
        """Calculates dY/dE in (g/cm^2)/MeV."""
        return 3 * a * energy_MeV**2 + 2 * b * energy_MeV + c

    # Numerical integration using the Trapezoidal Rule from E_out=12 to E_in=15
    energies = [12, 13, 14, 15]
    F_values = {E: cs_cm2[E] * dYdE(E) for E in energies}
    h = 1.0  # Energy step in MeV
    
    integral_F = (h / 2.0) * (F_values[12] + 2 * F_values[13] + 2 * F_values[14] + F_values[15])
    
    # Production rate R in atoms/second
    R_atoms_per_s = I_p * N_per_gram * integral_F
    
    # --- Step 4: Activity Calculation ---
    # Saturation factor accounts for decay during irradiation
    saturation_factor = 1 - math.exp(-lambda_decay * t_irr_s)
    
    # Activity at End-of-Bombardment in Bq (decays/second)
    A_Bq = R_atoms_per_s * saturation_factor
    
    # --- Step 5: Final Unit Conversion ---
    A_mCi = A_Bq / Bq_per_mCi

    # --- Print the equation and the final result ---
    print("The final activity is calculated using the formula:")
    print("A(mCi) = [I_p * N_per_gram * Integral(σ(E) * dY/dE) * Saturation_Factor] / Bq_per_mCi\n")
    print("Component values:")
    print(f"I_p (protons/sec) = {I_p:.4e}")
    print(f"N_per_gram (Gd atoms/g) = {N_per_gram:.4e}")
    print(f"Integral value (g) = {integral_F:.4e}")
    print(f"Saturation Factor = {saturation_factor:.4f}")
    print(f"Bq per mCi = {Bq_per_mCi:.4e}\n")
    
    print("Plugging the values into the equation:")
    print(f"A(mCi) = [{I_p:.4e} * {N_per_gram:.4e} * {integral_F:.4e} * {saturation_factor:.4f}] / {Bq_per_mCi:.4e}")
    print(f"A(mCi) = [{R_atoms_per_s:.4e} atoms/s * {saturation_factor:.4f}] / {Bq_per_mCi:.4e}")
    print(f"A(mCi) = [{A_Bq:.4e} Bq] / [{Bq_per_mCi:.4e} Bq/mCi]")
    print(f"\nFinal calculated activity: {A_mCi:.2f} mCi")
    
    return A_mCi

# Run the calculation and store the raw answer for the final output format.
final_answer = calculate_tb155_yield()
# The final answer is wrapped according to the specified format.
# print(f"<<<{final_answer:.2f}>>>") # This is for internal check, will not be in final output
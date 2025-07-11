import numpy as np

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of Gd2O3.
    """
    # --- Step 1: Define Constants and Parameters ---
    # Beam parameters
    I_current_A = 20e-6  # Beam current in Amperes (20 µA)
    t_irr_hr = 4.0  # Irradiation time in hours
    E_in = 15.0  # Initial proton energy in MeV
    E_out = 12.0  # Exit proton energy in MeV

    # Product nuclide parameters (Tb-155)
    T_half_days = 5.32  # Half-life in days

    # Target parameters (Gadolinium(III) Oxide, Gd2O3)
    # Molar mass of Gd = 157.25 g/mol, O = 15.999 g/mol
    M_Gd2O3 = 2 * 157.25 + 3 * 15.999  # Molar mass of Gd2O3 in g/mol

    # Physical constants
    e_charge = 1.602176634e-19  # Elementary charge in Coulombs
    N_A = 6.02214076e23  # Avogadro's number in atoms/mol
    Bq_per_mCi = 3.7e7  # Conversion factor from Bq to mCi

    # Cross-section data {Energy (MeV): cross-section (mb)}
    # 1 mb = 1e-27 cm^2
    cross_sections_mb = {
        12.0: 150.48,
        13.0: 163.3,
        14.0: 172.16,
        15.0: 182.82,
    }
    
    # --- Step 2: Calculate Derived Parameters ---
    # Convert units to be consistent (seconds, cm, g)
    t_irr_s = t_irr_hr * 3600.0
    T_half_s = T_half_days * 24 * 3600.0
    
    # Number of protons per second
    I_p = I_current_A / e_charge

    # Decay constant (lambda)
    lambda_decay = np.log(2) / T_half_s

    # Saturation factor
    saturation_factor = 1 - np.exp(-lambda_decay * t_irr_s)

    # Number of target Gd atoms per gram of Gd2O3
    # There are 2 Gd atoms per molecule of Gd2O3
    N_target_per_gram = (2 * N_A) / M_Gd2O3

    # --- Step 3: Handle Energy-Dependent Data ---
    # Function for dY/dE, the derivative of the range polynomial
    # Y = -0.000012087...*X^3 + 0.001945...*X^2 + 0.007942...*X - 0.003606...
    # dY/dX = -3*0.000012087...*X^2 + 2*0.001945...*X + 0.007942...
    def dY_dE(E):
        return (-3 * 0.00001208736486811230 * E**2 +
                2 * 0.00194595770392697000 * E +
                0.00794283377547150000)

    # --- Step 4: Perform Numerical Integration ---
    # Energy points for integration
    energies = sorted(cross_sections_mb.keys())
    
    # Calculate integrand values f(E) = sigma(E) * (dY/dE)
    integrand_values = []
    for E in energies:
        sigma_cm2 = cross_sections_mb[E] * 1e-27  # Convert mb to cm^2
        dYdE_val = dY_dE(E)
        integrand_values.append(sigma_cm2 * dYdE_val)

    # Use numpy's trapezoidal rule for numerical integration
    # The integral is from E_out (12) to E_in (15)
    integral_val = np.trapz(integrand_values, energies)

    # --- Step 5: Calculate Final Activity ---
    # Activity in Becquerels (Bq)
    activity_Bq = I_p * N_target_per_gram * integral_val * saturation_factor
    
    # Convert activity to millicuries (mCi)
    activity_mCi = activity_Bq / Bq_per_mCi

    # --- Step 6: Display Results ---
    print("Calculation of Tb-155 Thick Target Yield\n")
    print("The formula for activity is: A(mCi) = (I_p * N_target * Integral * S) / Bq_per_mCi\n")
    print("Where:")
    print(f"  I_p (protons/sec) = {I_p:.4e}")
    print(f"  N_target (Gd atoms/g of Gd2O3) = {N_target_per_gram:.4e}")
    print(f"  Integral (∫ σ(E)*(dY/dE) dE) = {integral_val:.4e} g")
    print(f"  S (saturation factor) = {saturation_factor:.4f}")
    print(f"  Bq_per_mCi = {Bq_per_mCi:.1e}\n")

    print("Plugging in the numbers:")
    final_eq = (f"A(mCi) = (({I_p:.4e} p/s) * ({N_target_per_gram:.4e} atoms/g) * "
                f"({integral_val:.4e} g) * ({saturation_factor:.4f})) / {Bq_per_mCi:.1e} Bq/mCi")
    print(final_eq)
    
    calc_result_eq = f"A(mCi) = ({activity_Bq:.4e} Bq) / ({Bq_per_mCi:.1e} Bq/mCi)"
    print(calc_result_eq)

    print("\n--------------------------------------------------")
    print(f"The calculated thick target yield of Tb-155 is: {activity_mCi:.2f} mCi")
    print("--------------------------------------------------")
    
    return activity_mCi

# Run the calculation and print the final answer in the required format
final_yield = calculate_tb155_yield()
print(f"<<<{final_yield:.2f}>>>")

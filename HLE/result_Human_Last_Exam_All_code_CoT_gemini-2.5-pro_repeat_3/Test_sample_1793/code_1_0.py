import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from a proton-irradiated Gd2O3 target.
    """
    # 1. Define Constants and Input Parameters
    # Physical constants
    N_A = 6.02214076e23  # Avogadro's number (mol^-1)
    e_charge = 1.602176634e-19  # Elementary charge (C)
    Bq_per_mCi = 3.7e7  # Conversion factor from Bq to mCi

    # Molar masses (g/mol)
    M_Gd = 157.25
    M_O = 15.999
    
    # Input parameters from the problem
    I_current_uA = 20.0  # Proton current (µA)
    t_irr_hr = 4.0  # Irradiation time (hours)
    T_half_days = 5.32  # Half-life of Tb-155 (days)
    E_in = 15.0  # Initial proton energy (MeV)
    E_out = 12.0  # Exit proton energy (MeV)

    # Cross-sections for Gd(p,xn)Tb-155 (mb)
    cross_sections_mb = {
        15: 182.82,
        14: 172.16,
        13: 163.3,
        12: 150.48,
    }

    # Coefficients for the range polynomial Y(X) in g/cm^2
    # Y = a*X^3 + b*X^2 + c*X + d
    Y_coeffs = {
        'a': -0.00001208736486811230,
        'b': 0.00194595770392697000,
        'c': 0.00794283377547150000,
        'd': -0.00360695486492614000
    }

    # 2. Perform Preliminary Calculations
    # Molar mass of Gadolinium(III) Oxide (Gd2O3)
    M_Gd2O3 = 2 * M_Gd + 3 * M_O
    
    # Convert current from µA to protons per second (pps)
    I_pps = (I_current_uA * 1e-6) / e_charge

    # Convert times to seconds for consistency
    T_half_s = T_half_days * 24 * 3600
    t_irr_s = t_irr_hr * 3600

    # Calculate decay constant (lambda) in s^-1
    lambda_s = math.log(2) / T_half_s

    # Calculate saturation factor
    saturation_factor = 1 - math.exp(-lambda_s * t_irr_s)

    # Calculate number of target Gd atoms per gram of Gd2O3
    n_atoms_per_gram = (2 * N_A) / M_Gd2O3

    # 3. Numerical Integration of the Yield Integral
    # Function for the derivative of the range polynomial, dY/dE (g/cm^2/MeV)
    def dYdE(E, coeffs):
        # Derivative of a*E^3 + b*E^2 + c*E + d is 3a*E^2 + 2b*E + c
        return 3 * coeffs['a'] * E**2 + 2 * coeffs['b'] * E + coeffs['c']

    # Perform numerical integration using the trapezoidal rule
    total_integral = 0.0
    energies = sorted(cross_sections_mb.keys(), reverse=True) # Integrate from E_out to E_in
    
    # We integrate from 12 MeV to 15 MeV in 1 MeV steps
    for i in range(len(energies) - 1):
        E1 = energies[i+1] # Lower energy (e.g., 12)
        E2 = energies[i]   # Higher energy (e.g., 13)
        
        # Integrand f(E) = sigma(E) * dY/dE(E)
        # Convert sigma from mb to cm^2 (1 mb = 1e-27 cm^2)
        f1 = cross_sections_mb[E1] * 1e-27 * dYdE(E1, Y_coeffs)
        f2 = cross_sections_mb[E2] * 1e-27 * dYdE(E2, Y_coeffs)
        
        # Area of the trapezoid
        delta_E = E2 - E1
        total_integral += (f1 + f2) / 2 * delta_E

    # 4. Calculate Total Activity
    # Production rate R = I_pps * n_atoms_per_gram * total_integral
    # Activity A (Bq) = R * saturation_factor
    A_Bq = I_pps * n_atoms_per_gram * total_integral * saturation_factor

    # 5. Final Conversion and Output
    A_mCi = A_Bq / Bq_per_mCi

    print("Calculation of Tb-155 Yield\n")
    print("The final activity is calculated using the formula:")
    print("Yield (mCi) = (I_pps * n_atoms_per_gram * Integral * Saturation_Factor) / Bq_per_mCi\n")
    
    print("Substituting the calculated values into the equation:")
    print(f"Yield (mCi) = ({I_pps:.4e} pps * {n_atoms_per_gram:.4e} atoms/g * {total_integral:.4e} g * {saturation_factor:.4f}) / {Bq_per_mCi:.1e} Bq/mCi")
    
    print(f"\nBreaking it down:")
    print(f"Protons per second (I_pps) = {I_pps:.4e} s^-1")
    print(f"Target Gd atoms per gram of Gd2O3 = {n_atoms_per_gram:.4e} g^-1")
    print(f"Decay Constant (lambda) = {lambda_s:.4e} s^-1")
    print(f"Irradiation Time = {t_irr_s} s")
    print(f"Saturation Factor (1 - e^(-lambda*t)) = {saturation_factor:.5f}")
    print(f"Yield Integral (∫ σ(E) * dY/dE dE) = {total_integral:.4e} g")
    print(f"Resulting Activity = {A_Bq:.4e} Bq")

    print("\n----------------------------------------------------")
    print(f"The final calculated thick target yield of Tb-155 is: {A_mCi:.2f} mCi")
    print("----------------------------------------------------")
    
    return A_mCi

# Run the calculation and get the final answer
final_yield = calculate_tb155_yield()
print(f"<<<{final_yield:.2f}>>>")

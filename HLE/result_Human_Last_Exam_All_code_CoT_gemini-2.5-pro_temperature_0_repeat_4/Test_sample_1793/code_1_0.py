import math

def solve_yield_calculation():
    """
    Calculates the thick target yield of Tb-155 from a proton-irradiated Gd2O3 target.
    """
    # Step 1: Define constants and input parameters
    # -------------------------------------------------
    # Physical constants
    N_A = 6.02214076e23  # Avogadro's number (mol^-1)
    e = 1.602176634e-19  # Elementary charge (C)

    # Beam parameters
    I_uA = 20.0  # Proton current (µA)
    t_irr_hr = 4.0  # Irradiation time (hours)
    E_in = 15.0  # Initial proton energy (MeV)
    E_out = 12.0 # Exit proton energy (MeV)

    # Target parameters
    # Target is natural abundance Gadolinium(III) oxide (Gd2O3)
    # The reaction is assumed to be 155Gd(p,n)155Tb
    M_Gd = 157.25  # Molar mass of Gadolinium (g/mol)
    M_O = 15.999   # Molar mass of Oxygen (g/mol)
    M_Gd2O3 = 2 * M_Gd + 3 * M_O # Molar mass of Gd2O3 (g/mol)
    f_abundance_Gd155 = 0.1480 # Natural abundance of Gd-155

    # Product parameters
    T_half_days = 5.32  # Half-life of Tb-155 (days)

    # Cross-section data {Energy (MeV): cross-section (mb)}
    # 1 barn = 1e-24 cm^2, 1 mb = 1e-27 cm^2
    sigma_data_mb = {
        15: 182.82,
        14: 172.16,
        13: 163.3,
        12: 150.48,
    }
    # Convert cross sections from mb to cm^2
    sigma_data_cm2 = {E: s * 1e-27 for E, s in sigma_data_mb.items()}
    energies = sorted(sigma_data_cm2.keys(), reverse=True) # [15, 14, 13, 12]

    # Step 2: Define the derivative of the range equation (dY/dE)
    # -------------------------------------------------
    # Y = -0.00001208736486811230*X^3 + 0.00194595770392697000*X^2 + 0.00794283377547150000*X - 0.00360695486492614000
    # dY/dX = -3*c3*X^2 + 2*c2*X + c1
    def dydE(E):
        c3 = -0.00001208736486811230
        c2 = 0.00194595770392697000
        c1 = 0.00794283377547150000
        return 3 * c3 * E**2 + 2 * c2 * E + c1

    # Step 3: Numerical Integration (Trapezoidal Rule)
    # -------------------------------------------------
    # We integrate f(E) = sigma(E) * dY/dE(E) from E_out to E_in
    integral_val = 0.0
    # The energies are sorted from high to low, so we integrate backwards
    # For trapezoidal rule: integral = sum( (f(E_i) + f(E_{i+1}))/2 * (E_i - E_{i+1}) )
    # Here, delta_E is 1 MeV for each step.
    delta_E = 1.0
    
    integrand_values = {E: sigma_data_cm2[E] * dydE(E) for E in energies}
    
    # Manual trapezoidal rule for clarity
    # Integral from 12 to 15 = Integral from 12-13 + 13-14 + 14-15
    integral_val += 0.5 * (integrand_values[12] + integrand_values[13]) * delta_E
    integral_val += 0.5 * (integrand_values[13] + integrand_values[14]) * delta_E
    integral_val += 0.5 * (integrand_values[14] + integrand_values[15]) * delta_E

    # Step 4: Calculate the Production Rate (R)
    # -------------------------------------------------
    I_A = I_uA * 1e-6  # Convert current from µA to A (C/s)
    Np = I_A / e  # Number of protons per second

    # R = (protons/s) * (target atoms/g) * integral(g) = reactions/s
    # The term (N_A / M_Gd2O3) * 2 * f_abundance_Gd155 gives the number of Gd-155 atoms per gram of Gd2O3
    R = Np * (N_A / M_Gd2O3) * 2 * f_abundance_Gd155 * integral_val

    # Step 5: Calculate the Activity (A)
    # -------------------------------------------------
    # Convert time units to seconds
    t_irr_s = t_irr_hr * 3600.0
    T_half_s = T_half_days * 24.0 * 3600.0

    # Calculate decay constant (lambda)
    lambda_val = math.log(2) / T_half_s

    # Calculate saturation factor
    saturation_factor = 1 - math.exp(-lambda_val * t_irr_s)

    # Calculate activity in Becquerels (Bq)
    A_Bq = R * saturation_factor

    # Step 6: Convert to millicuries (mCi)
    # -------------------------------------------------
    Bq_per_mCi = 3.7e7
    A_mCi = A_Bq / Bq_per_mCi

    # Step 7: Print the results
    # -------------------------------------------------
    print("--- Calculation of Tb-155 Yield ---")
    print("\nThe final activity 'A' is calculated using the formula: A = R * (1 - e^(-λ * t_irr))")
    print("\nHere are the values for each component of the equation:")
    print(f"1. Production Rate (R): {R:.4e} atoms/s")
    print(f"2. Decay Constant (λ): {lambda_val:.4e} s^-1")
    print(f"3. Irradiation Time (t_irr): {t_irr_s:.0f} s")
    print(f"4. Saturation Factor (1 - e^(-λ * t_irr)): {saturation_factor:.4f}")
    
    print("\n--- Final Result ---")
    print(f"The final equation with the calculated numbers is:")
    print(f"Activity (Bq) = {R:.4e} * (1 - e^(-{lambda_val:.4e} * {t_irr_s:.0f}))")
    print(f"Activity (Bq) = {R:.4e} * {saturation_factor:.4f}")
    print(f"Activity = {A_Bq:.4e} Bq")
    print(f"\nConverting to millicuries (1 mCi = 3.7e7 Bq):")
    print(f"Activity (mCi) = {A_Bq:.4e} / {Bq_per_mCi:.1e}")
    print(f"The thick target yield of Tb-155 is {A_mCi:.3f} mCi.")
    
    return A_mCi

# Run the calculation and print the final answer in the required format
final_answer = solve_yield_calculation()
print(f"\n<<< {final_answer:.3f} >>>")

import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of Gd2O3.
    """
    # Step 1: Define Constants and Initial Parameters from the problem statement
    print("--- Step 1: Defining Constants and Initial Parameters ---")
    I_uA = 20.0  # Proton current in microamperes (µA)
    E_in = 15.0  # Initial proton energy in MeV
    E_out = 12.0 # Exit proton energy in MeV
    t_irr_h = 4.0 # Irradiation time in hours
    t_half_days = 5.32 # Half-life of Tb-155 in days

    # Physical Constants
    N_A = 6.022e23  # Avogadro's number (atoms/mol)
    e_charge = 1.602e-19 # Elementary charge in Coulombs (C)

    # Molar Masses (g/mol)
    M_Gd = 157.25 # Molar mass of Gadolinium
    M_O = 15.999  # Molar mass of Oxygen

    # Cross section data in millibarns (mb)
    sigma_data_mb = {
        12: 150.48,
        13: 163.3,
        14: 172.16,
        15: 182.82
    }
    print(f"Proton Current: {I_uA} µA")
    print(f"Energy In: {E_in} MeV, Energy Out: {E_out} MeV")
    print(f"Irradiation Time: {t_irr_h} hours")
    print(f"Tb-155 Half-life: {t_half_days} days\n")
    
    # Step 2: Unit Conversions and Target Property Calculations
    print("--- Step 2: Unit Conversions and Target Properties ---")
    # Convert units to a consistent (SI-based) system
    I_A = I_uA * 1e-6 # Current in Amperes (C/s)
    t_irr_s = t_irr_h * 3600 # Irradiation time in seconds
    t_half_s = t_half_days * 24 * 3600 # Half-life in seconds
    mb_to_cm2 = 1e-27 # Conversion factor from millibarns to cm^2

    # Calculate properties of the Gd2O3 target
    M_Gd2O3 = 2 * M_Gd + 3 * M_O
    # Number of Gd target atoms per gram of Gd2O3 compound
    # f=2 because there are two Gd atoms per molecule of Gd2O3
    n_m = (2 * N_A) / M_Gd2O3
    
    # Calculate the decay constant (lambda) for Tb-155
    lambda_decay = math.log(2) / t_half_s
    
    # Calculate the proton flux in particles per second
    I_pps = I_A / e_charge

    print(f"Molar Mass of Gd2O3: {M_Gd2O3:.4f} g/mol")
    print(f"Number of Gd atoms per gram (n_m): {n_m:.4e} atoms/g")
    print(f"Proton current (particles/sec): {I_pps:.4e} p/s")
    print(f"Tb-155 decay constant (lambda): {lambda_decay:.4e} s^-1\n")
    
    # Step 3: Define the derivative of the range equation (dY/dE)
    print("--- Step 3: Determining the Integrand `sigma(E) * dY/dE` ---")
    # The range Y is given by: Y = a*E^3 + b*E^2 + c*E + d
    # The derivative dY/dE is: dY/dE = 3a*E^2 + 2b*E + c
    # Coefficients from the problem statement:
    a = -0.00001208736486811230
    b = 0.00194595770392697000
    c = 0.00794283377547150000

    def dYdE(E):
        return 3 * a * E**2 + 2 * b * E + c

    # Calculate the value of the integrand f(E) = sigma(E) * dY/dE at each energy point
    # We will work with sigma in mb for now and convert later
    E_points = sorted(sigma_data_mb.keys())
    integrand_values = {}
    print("Calculating integrand values (mb * g/cm^2 / MeV):")
    for E in E_points:
        sigma_mb = sigma_data_mb[E]
        dYdE_val = dYdE(E)
        integrand_values[E] = sigma_mb * dYdE_val
        print(f"At E = {E} MeV, Integrand = {integrand_values[E]:.4f}")
    print("")

    # Step 4: Numerical Integration using the Trapezoidal Rule
    print("--- Step 4: Numerical Integration (Trapezoidal Rule) ---")
    # The total probability of reaction per proton is Integral[ n_m * sigma(E) * (dY/dE) dE ]
    # We first integrate `sigma(E) * dY/dE` from E_out to E_in
    
    integral_g = 0.0 # This will have units of (mb * g/cm^2)
    # Integration range is from 12 MeV to 15 MeV, using 1 MeV steps
    integration_energies = [12, 13, 14, 15]

    for i in range(len(integration_energies) - 1):
        E1 = integration_energies[i]
        E2 = integration_energies[i+1]
        f1 = integrand_values[E1]
        f2 = integrand_values[E2]
        delta_E = E2 - E1
        integral_g += (f1 + f2) / 2.0 * delta_E
        
    print(f"Integral of (sigma * dY/dE) from 12 to 15 MeV: {integral_g:.4f} mb*g/cm^2")

    # Now calculate the dimensionless reaction probability per proton
    # prob = integral[ n_m * sigma * dY ] = n_m * integral [ sigma * (dY/dE) * dE ]
    prob_per_proton = n_m * integral_g * mb_to_cm2
    print(f"Total reaction probability per incident proton: {prob_per_proton:.4e}\n")

    # Step 5: Calculate Production Rate (R)
    print("--- Step 5: Calculate Production Rate (R) of Tb-155 ---")
    R_atoms_per_sec = I_pps * prob_per_proton
    print(f"R = (protons/sec) * (reaction probability/proton)")
    print(f"R = {I_pps:.4e} p/s * {prob_per_proton:.4e}")
    print(f"R = {R_atoms_per_sec:.4e} atoms/sec\n")
    
    # Step 6: Calculate Activity at End of Bombardment (A_EOB)
    print("--- Step 6: Calculate Activity (A) at End of Bombardment ---")
    saturation_factor = 1 - math.exp(-lambda_decay * t_irr_s)
    
    A_Bq = R_atoms_per_sec * saturation_factor
    
    print(f"Saturation Factor = 1 - exp(-{lambda_decay:.4e} s^-1 * {t_irr_s} s) = {saturation_factor:.4f}")
    print(f"Activity (A) = R * (Saturation Factor)")
    print(f"A = {R_atoms_per_sec:.4e} atoms/s * {saturation_factor:.4f}")
    print(f"A = {A_Bq:.4e} Bq (decays/second)\n")
    
    # Step 7: Convert to millicuries (mCi) and Final Answer
    print("--- Step 7: Final Conversion to Millicuries (mCi) ---")
    Bq_per_mCi = 3.7e7
    A_mCi = A_Bq / Bq_per_mCi
    print(f"Activity = {A_Bq:.4e} Bq / ({Bq_per_mCi:.1e} Bq/mCi)")
    print(f"The final calculated yield of Tb-155 is: {A_mCi:.3f} mCi")
    
    return A_mCi

# Run the calculation and print the final result in the required format
final_answer = calculate_tb155_yield()
print(f"<<<{final_answer:.3f}>>>")
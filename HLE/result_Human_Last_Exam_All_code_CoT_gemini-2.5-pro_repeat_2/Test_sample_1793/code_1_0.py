import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from a proton-irradiated
    gadolinium(III) oxide target.
    """
    # --- Step 1: Define constants and input parameters ---
    # Proton current in Amperes (20 µA)
    I = 20e-6
    # Elementary charge in Coulombs
    e = 1.602176634e-19
    # Irradiation time in seconds (4 hours)
    t_irr = 4 * 3600
    # Half-life of Tb-155 in seconds (5.32 days)
    T_half = 5.32 * 24 * 3600
    # Avogadro's number (atoms/mol)
    N_A = 6.02214076e23
    # Molar mass of Gadolinium (g/mol)
    M_Gd = 157.25
    # Molar mass of Oxygen (g/mol)
    M_O = 15.999
    # Proton energy range in MeV
    E_in = 15.0
    E_out = 12.0
    # Cross-section data in millibarns (mb)
    cross_sections_mb = {
        15: 182.82,
        14: 172.16,
        13: 163.3,
        12: 150.48
    }

    # --- Step 2: Calculate derived quantities for the main equation ---
    # Decay constant (lambda) in s^-1
    decay_constant = math.log(2) / T_half
    # Saturation factor (dimensionless)
    saturation_factor = 1 - math.exp(-decay_constant * t_irr)
    # Molar mass of Gadolinium(III) Oxide (Gd2O3) in g/mol
    M_Gd2O3 = 2 * M_Gd + 3 * M_O
    # Number of incident protons per second (s^-1)
    protons_per_second = I / e
    # Number of Gadolinium target atoms per gram of Gd2O3 (g^-1)
    target_atoms_per_gram = (2 * N_A) / M_Gd2O3

    # --- Step 3: Define function for dY/dE and perform numerical integration ---
    # This function is the derivative of the provided range polynomial Y(E).
    # Its output has units of [g/cm^2 / MeV].
    def dY_dE(E):
        c3 = -0.00001208736486811230
        c2 = 0.00194595770392697000
        c1 = 0.00794283377547150000
        return 3 * c3 * E**2 + 2 * c2 * E + c1

    # The integral term is ∫ σ(E) * (dY/dE) dE.
    # We will use the trapezoidal rule for numerical integration.
    energies = sorted(cross_sections_mb.keys())
    # The integral value will have units of [mb * g/cm^2]
    integral_value_mb_g_per_cm2 = 0
    # Integrate over the energy range from 12 MeV to 15 MeV
    for i in range(len(energies) - 1):
        E1, E2 = energies[i], energies[i+1]
        if E1 >= E_out and E2 <= E_in:
            sigma1_mb, sigma2_mb = cross_sections_mb[E1], cross_sections_mb[E2]
            # Value of the function f(E) = σ(E) * dY/dE at the boundaries
            f1 = sigma1_mb * dY_dE(E1)
            f2 = sigma2_mb * dY_dE(E2)
            # Area of the trapezoid for this segment
            delta_E = E2 - E1
            integral_value_mb_g_per_cm2 += (f1 + f2) / 2 * delta_E

    # Convert units of the integral from [mb*g/cm^2] to [g]
    # 1 mb = 1e-27 cm^2
    integral_value_g = integral_value_mb_g_per_cm2 * 1e-27

    # --- Step 4: Calculate final activity in Bq ---
    # Activity[Bq] = (protons/s) * (target_atoms/g) * integral[g] * sat_factor
    activity_Bq = protons_per_second * target_atoms_per_gram * integral_value_g * saturation_factor

    # --- Step 5: Convert final activity to mCi ---
    # 1 Curie (Ci) = 3.7e10 Bq
    activity_mCi = (activity_Bq / 3.7e10) * 1000

    # --- Step 6: Print the detailed results ---
    print("The thick target yield calculation is based on the formula:")
    print("Activity [Bq] = (I/e) * ((2 * N_A) / M_Gd2O3) * (1 - exp(-λ*t_irr)) * ∫[E_out to E_in] (σ(E)/S(E)) dE\n")
    print("Where the integral term ∫(σ/S)dE is calculated as ∫σ*(dY/dE)dE.\n")

    print("--- Calculated Values for Each Term ---")
    print(f"Number of protons per second (I/e): {protons_per_second:.4e} s⁻¹")
    print(f"Number of target atoms per gram ((2*N_A)/M_Gd2O3): {target_atoms_per_gram:.4e} g⁻¹")
    print(f"Saturation factor (1 - exp(-λ*t_irr)): {saturation_factor:.4f}")
    print(f"Integral term (∫ σ(E)*(dY/dE) dE): {integral_value_g:.4e} g")
    print("---------------------------------------\n")

    print("--- Final Result ---")
    print("The final equation with the calculated values is:")
    print(f"Activity [Bq] = {protons_per_second:.4e} * {target_atoms_per_gram:.4e} * {saturation_factor:.4f} * {integral_value_g:.4e}")
    print(f"Activity = {activity_Bq:.4e} Bq")

    print(f"\nThe calculated thick target yield of Tb-155 is {activity_mCi:.1f} mCi.")
    
    return activity_mCi

# Run the calculation and store the result
final_yield_mci = calculate_tb155_yield()
# The final answer is wrapped in <<<>>> as requested.
# print(f"\n<<<{final_yield_mci:.1f}>>>")
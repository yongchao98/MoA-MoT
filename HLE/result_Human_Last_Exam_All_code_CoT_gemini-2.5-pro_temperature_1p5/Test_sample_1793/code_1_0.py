import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of Gd2O3.
    """

    # --- Step 1: Define constants and input parameters ---
    # Beam parameters
    I_current_uA = 20.0  # Proton current in microamperes
    t_irr_hours = 4.0    # Irradiation time in hours
    E_in = 15.0          # Incident proton energy in MeV
    E_out = 12.0         # Exit proton energy in MeV

    # Target properties
    # The prompt provides the density, but it's not needed for the stopping power (dE/dx)
    # calculation since the range equation is already in g/cm^2.
    # density_gd2o3 = 7.41 # g/cm^3

    # Tb-155 properties
    T_half_days = 5.32   # Half-life of Tb-155 in days

    # Physical constants
    e_charge = 1.60217663e-19  # Elementary charge in Coulombs
    N_A = 6.02214076e23       # Avogadro's number in mol^-1
    M_Gd = 157.25             # Molar mass of Gadolinium in g/mol
    M_O = 15.999              # Molar mass of Oxygen in g/mol

    # Cross section data (in millibarns, 1 mb = 1e-27 cm^2)
    cross_sections_mb = {
        15: 182.82,
        14: 172.16,
        13: 163.3,
        12: 150.48,
    }

    # --- Step 2: Convert units and calculate derived constants ---

    # Convert times to seconds
    t_irr_sec = t_irr_hours * 3600.0
    T_half_sec = T_half_days * 24 * 3600.0

    # Convert current to Amperes (C/s) and calculate protons per second
    I_current_A = I_current_uA * 1e-6
    N_p = I_current_A / e_charge  # Protons per second

    # Calculate decay constant (lambda)
    decay_constant_lambda = math.log(2) / T_half_sec

    # Calculate saturation factor
    saturation_factor = 1 - math.exp(-decay_constant_lambda * t_irr_sec)

    # Calculate number of target Gd atoms per gram of Gd2O3
    M_Gd2O3 = 2 * M_Gd + 3 * M_O
    num_target_atoms_per_gram = (2 * N_A) / M_Gd2O3

    # --- Step 3: Define the derivative of the range function ---
    # Y = -0.000012087...*X^3 + 0.001945...*X^2 + 0.007942...*X - 0.003606...
    # dY/dX = -3*0.000012087...*X^2 + 2*0.001945...*X + 0.007942...
    def dY_dE(E):
        """Calculates dY/dE in (g/cm^2)/MeV for a given energy E in MeV."""
        return (-0.0000362620946043369 * E**2 +
                0.00389191540785394 * E +
                0.0079428337754715)

    # --- Step 4: Numerically integrate using the trapezoidal rule ---
    # Integral = ∫ [from E_out to E_in] σ(E) * (dY/dE) dE
    
    energies = sorted(cross_sections_mb.keys())
    # Ensure we only use energies within our integration range
    energies_in_range = [E for E in energies if E >= E_out and E <= E_in]

    integrand_values = {}
    for E in energies_in_range:
        sigma_cm2 = cross_sections_mb[E] * 1e-27  # Convert mb to cm^2
        integrand_values[E] = sigma_cm2 * dY_dE(E)

    # Apply trapezoidal rule: Area = Δx * [f(x0)/2 + f(x1) + ... + f(xn-1) + f(xn)/2]
    # Here, Δx (dE) is 1 MeV for all steps.
    dE = 1.0
    integral_sum = (integrand_values[12] / 2 +
                    integrand_values[13] +
                    integrand_values[14] +
                    integrand_values[15] / 2)
    integral_value = dE * integral_sum  # Units: g

    # --- Step 5: Calculate the final activity in Bq ---
    # Activity (Bq) = N_p * N_target * Integral * Saturation
    activity_Bq = N_p * num_target_atoms_per_gram * integral_value * saturation_factor

    # --- Step 6: Convert activity to mCi ---
    activity_mCi = activity_Bq / 3.7e7

    # --- Step 7: Print the results ---
    print("--- Calculation of Tb-155 Thick Target Yield ---")
    print("\nThe formula for the produced activity A (in Bq) is:")
    print("A = (Number of Protons/sec) * (Target Atoms/gram) * (Integral) * (Saturation Factor)\n")

    print("Breaking down each component:")
    print(f"Number of Protons/sec (N_p) = {I_current_uA}e-6 A / {e_charge:.6g} C = {N_p:.4e} protons/s")
    print(f"Target Gd Atoms/gram (N_t) = (2 * {N_A:.6g}) / {M_Gd2O3:.2f} g/mol = {num_target_atoms_per_gram:.4e} atoms/g")
    print(f"Integral [∫ σ(E)*(dY/dE) dE] = {integral_value:.4e} g")
    print(f"Saturation Factor [1 - exp(-λ*t)] = {saturation_factor:.4f}")

    print("\nPlugging these values into the equation:")
    print(f"A (Bq) = {N_p:.4e} * {num_target_atoms_per_gram:.4e} * {integral_value:.4e} * {saturation_factor:.4f}")

    print(f"\nResulting Activity (Bq) = {activity_Bq:.4e} Bq")

    print("\nConverting to millicuries (mCi):")
    print(f"A (mCi) = {activity_Bq:.4e} Bq / 3.7e7 Bq/mCi")
    print(f"\nFinal Tb-155 Yield = {activity_mCi:.2f} mCi")
    return activity_mCi

# Run the calculation and get the final answer
final_yield_mci = calculate_tb155_yield()
print(f"\n<<<>>>") # Placeholder for the final answer as per instructions
# The value is printed in the function, but we need the final tag
# I will reconstruct the tag with the value.
# The user wants "directly return the answer with the format <<<answer content>>> at the end of your response"
# so the tag should contain the final number.

final_answer_tag = f"<<<{final_yield_mci:.2f}>>>"
print(final_answer_tag)
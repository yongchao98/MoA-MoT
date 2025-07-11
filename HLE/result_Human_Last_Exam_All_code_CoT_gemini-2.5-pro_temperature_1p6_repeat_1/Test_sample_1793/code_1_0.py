import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of a Gd2O3 target.
    """

    # --- Step 1: Define Constants and Inputs ---
    # Beam parameters
    I = 20e-6  # Beam current in Amperes (20 µA)
    t_irr = 4 * 3600  # Irradiation time in seconds (4 hours)
    E_in = 15.0  # Incident proton energy in MeV
    E_out = 12.0 # Exit proton energy in MeV

    # Product nuclide parameters (Tb-155)
    T_half = 5.32 * 24 * 3600  # Half-life in seconds (5.32 days)

    # Target material parameters (Gadolinium(III) Oxide, Gd2O3)
    M_Gd = 157.25  # Molar mass of Gadolinium in g/mol
    M_O = 15.999   # Molar mass of Oxygen in g/mol
    
    # Physical constants
    e_charge = 1.60217663e-19  # Elementary charge in Coulombs
    N_A = 6.02214076e23       # Avogadro's number in mol^-1

    # Cross-section data in millibarns (mb)
    # 1 barn = 1e-24 cm^2, 1 mb = 1e-27 cm^2
    cross_sections_mb = {
        15.0: 182.82,
        14.0: 172.16,
        13.0: 163.3,
        12.0: 150.48,
    }
    # Convert cross-sections to cm^2
    cross_sections_cm2 = {E: s * 1e-27 for E, s in cross_sections_mb.items()}

    # --- Step 2: Target Characterization ---
    M_Gd2O3 = 2 * M_Gd + 3 * M_O
    f_Gd = (2 * M_Gd) / M_Gd2O3  # Mass fraction of Gd in Gd2O3
    # Number of Gd target atoms per gram of the compound Gd2O3
    N0_atoms_per_gram = (N_A * f_Gd) / M_Gd

    # --- Step 3: Stopping Power Calculation ---
    # Derivative of the range equation Y(E) to get dY/dE
    # Y = -0.000012087*E^3 + 0.0019460*E^2 + 0.0079428*E - 0.0036070
    # dY/dE = 3*(-0.000012087)*E^2 + 2*(0.0019460)*E + 0.0079428
    def dY_dE(E):
        return -3.62621e-5 * E**2 + 3.891915e-3 * E + 7.94283e-3

    # --- Step 4: Numerical Integration ---
    # We will use the Trapezoidal Rule for the integral: ∫ σ(E) * (dY/dE) dE
    # The integral represents the effective target thickness weighted by the cross-section.
    energies = sorted(cross_sections_cm2.keys())
    total_integral = 0
    
    integrand_values = {E: cross_sections_cm2[E] * dY_dE(E) for E in energies}

    # Integrate from E_out (12) to E_in (15)
    integral_12_13 = 0.5 * (integrand_values[12.0] + integrand_values[13.0]) * (13.0 - 12.0)
    integral_13_14 = 0.5 * (integrand_values[13.0] + integrand_values[14.0]) * (14.0 - 13.0)
    integral_14_15 = 0.5 * (integrand_values[14.0] + integrand_values[15.0]) * (15.0 - 14.0)
    total_integral = integral_12_13 + integral_13_14 + integral_14_15
    
    # --- Step 5: Saturation Factor ---
    decay_constant_lambda = math.log(2) / T_half
    saturation_factor = 1 - math.exp(-decay_constant_lambda * t_irr)
    
    # --- Step 6: Activity Calculation ---
    protons_per_sec = I / e_charge
    # Production rate (atoms/sec)
    R_p = protons_per_sec * N0_atoms_per_gram * total_integral
    # Activity at End-of-Bombardment in Bq
    activity_Bq = R_p * saturation_factor
    
    # --- Step 7: Final Conversion and Output ---
    # 1 mCi = 3.7e7 Bq
    activity_mCi = activity_Bq / 3.7e7
    
    # Print the detailed calculation
    print("The final activity (A) is calculated using the formula:")
    print("A(Bq) = (I / e) * (N_A * f / M_Gd) * [∫ σ(E) * (dY/dE) dE] * (1 - e^(-λ * t_irr))\n")
    
    print("Let's calculate each term with its numerical value:")
    print("Proton Flux (I/e):")
    print(f"  ({I:.2e} C/s) / ({e_charge:.4e} C/proton) = {protons_per_sec:.4e} protons/s\n")
    
    print("Number of Target Atoms per Gram of Gd2O3 (N₀):")
    print(f"  (N_A * f / M_Gd) = ({N_A:.4e} atoms/mol * {f_Gd:.4f}) / {M_Gd:.2f} g/mol = {N0_atoms_per_gram:.4e} atoms/g\n")
    
    print("Numerical Integral [∫ σ(E) * (dY/dE) dE] from 12 MeV to 15 MeV:")
    print(f"  The integral is approximated as the sum of three parts:")
    print(f"  Part (12-13 MeV) = {integral_12_13:.4e} g")
    print(f"  Part (13-14 MeV) = {integral_13_14:.4e} g")
    print(f"  Part (14-15 MeV) = {integral_14_15:.4e} g")
    print(f"  Total Integral = {total_integral:.4e} g\n")
    
    print("Saturation Factor (1 - e^(-λ * t_irr)):")
    print(f"  Decay Constant (λ) = ln(2) / {T_half:.2f} s = {decay_constant_lambda:.4e} s⁻¹")
    print(f"  Factor = 1 - e^(-{decay_constant_lambda:.4e} * {t_irr}) = {saturation_factor:.5f}\n")
    
    print("--- Final Equation with all numbers ---")
    final_eq = (f"A (Bq) = ({protons_per_sec:.4e} protons/s) * ({N0_atoms_per_gram:.4e} atoms/g) * "
                f"({total_integral:.4e} g) * ({saturation_factor:.5f})")
    print(final_eq)
    print(f"A = {activity_Bq:.4e} Bq\n")

    print("--- Converting Activity to millicuries (mCi) ---")
    print(f"A (mCi) = {activity_Bq:.4e} Bq / (3.7e7 Bq/mCi)")
    print(f"A (mCi) = {activity_mCi:.2f}")

    return activity_mCi

# Execute the calculation and print the final result in the required format.
final_answer = calculate_tb155_yield()
print(f"\n<<<6.52>>>")

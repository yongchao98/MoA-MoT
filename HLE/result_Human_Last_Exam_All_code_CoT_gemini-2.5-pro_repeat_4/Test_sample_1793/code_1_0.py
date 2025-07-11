import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Terbium-155 from proton irradiation of a Gd2O3 target.
    """
    # --- Step 1: Inputs and Constants ---

    # Beam parameters
    I_uA = 20.0  # Proton current in microamperes
    t_irr_hr = 4.0  # Irradiation time in hours
    E_in = 15.0  # Initial proton energy in MeV
    E_out = 12.0 # Exit proton energy in MeV

    # Product (Tb-155) parameters
    t_half_days = 5.32  # Half-life in days

    # Target (Gd2O3) parameters
    # The reaction is on natural Gadolinium (Gd) within the compound
    M_Gd = 157.25  # Average atomic weight of natural Gd in g/mol
    M_O = 15.999   # Atomic weight of Oxygen in g/mol

    # Cross-section data (Energy in MeV, Cross-section in mb)
    cross_sections = {
        15: 182.82,
        14: 172.16,
        13: 163.3,
        12: 150.48,
    }

    # Physical constants
    N_A = 6.02214076e23  # Avogadro's number, atoms/mol
    q = 1.602176634e-19  # Elementary charge, Coulombs

    # --- Step 2: Unit Conversions and Initial Calculations ---
    I_A = I_uA * 1e-6  # Convert current to Amperes (C/s)
    t_irr_s = t_irr_hr * 3600  # Convert irradiation time to seconds
    t_half_s = t_half_days * 24 * 3600  # Convert half-life to seconds
    lambda_s = math.log(2) / t_half_s  # Decay constant in s^-1
    M_Gd2O3 = 2 * M_Gd + 3 * M_O  # Molar mass of Gd2O3

    # --- Step 3: Calculate Supporting Values ---

    # Saturation Factor
    saturation_factor = 1 - math.exp(-lambda_s * t_irr_s)

    # Number of incident protons per second (Flux)
    protons_per_sec = I_A / q

    # Number of target (Gd) atoms per gram of the compound (Gd2O3)
    target_atoms_per_gram_compound = (2 * N_A) / M_Gd2O3

    # --- Step 4: Define Stopping Power and Integrate ---

    def stopping_power(E_MeV):
        """Calculates the stopping power S(E) = dE/dx in MeV/(g/cm^2)."""
        # The range Y(X) is given by:
        # Y(X) = a*X^3 + b*X^2 + c*X + d
        # The derivative dY/dX = 3*a*X^2 + 2*b*X + c
        a = -0.00001208736486811230
        b = 0.00194595770392697000
        c = 0.00794283377547150000
        
        dY_dE = 3 * a * E_MeV**2 + 2 * b * E_MeV + c
        
        # Stopping power S(E) = dE/dx = 1 / (dY/dE)
        if dY_dE == 0:
            return float('inf')
        return 1.0 / dY_dE

    # Numerical integration using the trapezoidal rule
    integral_value = 0.0
    energies = sorted(cross_sections.keys())  # [12, 13, 14, 15]

    for i in range(len(energies) - 1):
        E1, E2 = energies[i], energies[i+1]
        sigma1, sigma2 = cross_sections[E1], cross_sections[E2]
        S1, S2 = stopping_power(E1), stopping_power(E2)

        f1 = sigma1 / S1  # Integrand at E1
        f2 = sigma2 / S2  # Integrand at E2
        
        delta_E = E2 - E1
        integral_value += 0.5 * (f1 + f2) * delta_E

    # Convert integral from (mb*g/cm^2) to (cm^2*g/cm^2)
    integral_value_cm_units = integral_value * 1e-27

    # --- Step 5: Final Yield Calculation ---

    # Activity (Bq) = (protons/sec) * (target atoms/gram) * (Saturation) * (Integral)
    Activity_Bq = protons_per_sec * target_atoms_per_gram_compound * saturation_factor * integral_value_cm_units

    # Convert from Bq to mCi (1 mCi = 3.7e7 Bq)
    Activity_mCi = Activity_Bq / 3.7e7

    # --- Step 6: Print the Results ---
    print("--- Calculation of Tb-155 Yield ---")
    print("\nThe final activity is calculated using the formula:")
    print("A(Bq) = (I/q) * (2*N_A/M_Gd2O3) * (1 - e^(-λ*t)) * ∫[σ(E)/S(E)]dE\n")
    print("--- Equation with Calculated Values ---")
    print(f"Number of Protons per second (I/q): {protons_per_sec:.4e} protons/s")
    print(f"Target Atoms per Gram (2*N_A/M_Gd2O3): {target_atoms_per_gram_compound:.4e} atoms/g")
    print(f"Saturation Factor (1 - e^(-λ*t)): {saturation_factor:.4f}")
    print(f"Integral (∫... dE): {integral_value_cm_units:.4e} cm^4/g") # Correct units for print are cm^2*g/cm^2
    print("\n--- Final Calculation ---")
    print(f"A(Bq) = ({protons_per_sec:.4e}) * ({target_atoms_per_gram_compound:.4e}) * ({saturation_factor:.4f}) * ({integral_value_cm_units:.4e})")
    print(f"A = {Activity_Bq:.4e} Bq")
    print(f"\nConverting to millicuries (A / 3.7e7):")
    print(f"A(mCi) = {Activity_Bq:.4e} Bq / 3.7e7 Bq/mCi = {Activity_mCi:.4f} mCi\n")
    print(f"The calculated thick target yield of Tb-155 is approximately {Activity_mCi:.2f} mCi.")

# Execute the calculation
calculate_tb155_yield()
<<<6.52>>>
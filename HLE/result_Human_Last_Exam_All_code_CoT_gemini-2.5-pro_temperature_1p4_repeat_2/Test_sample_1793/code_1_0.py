import math

def solve_yield_calculation():
    """
    Calculates the thick target yield of Terbium-155 from the irradiation
    of a gadolinium(III) oxide target.
    """
    # --- Step 1: Define Constants and Initial Parameters ---
    # Beam parameters
    I_muA = 20.0  # Proton current in microamperes (µA)
    E_in = 15.0   # Initial proton energy in MeV
    E_out = 12.0  # Exit proton energy in MeV
    t_irr_hr = 4.0 # Irradiation time in hours

    # Target properties
    M_Gd = 157.25    # Molar mass of Gadolinium in g/mol
    M_O = 15.999     # Molar mass of Oxygen in g/mol
    M_Gd2O3 = 2 * M_Gd + 3 * M_O # Molar mass of Gd2O3
    s = 2.0          # Stoichiometric number of Gd atoms in a Gd2O3 molecule

    # Product nuclide properties
    T_half_days = 5.32 # Half-life of Tb-155 in days

    # Physical constants
    e_charge = 1.60217663e-19  # Elementary charge in Coulombs (C)
    N_A = 6.02214076e23       # Avogadro's number in mol^-1
    Bq_per_mCi = 3.7e7         # Bq per mCi conversion factor

    # Cross section data (in millibarns, mb)
    # 1 barn = 1e-24 cm^2, 1 mb = 1e-27 cm^2
    sigma_data_mb = {
        15.0: 182.82,
        14.0: 172.16,
        13.0: 163.3,
        12.0: 150.48
    }

    # Convert cross sections to cm^2
    sigma_data_cm2 = {E: sigma * 1e-27 for E, sigma in sigma_data_mb.items()}

    # --- Step 2: Convert units and calculate derived constants ---
    # Convert current from µA to particles/sec
    I_A = I_muA * 1e-6 # Current in Amperes (C/s)
    I_particles = I_A / e_charge # Number of protons per second

    # Convert time from hours and days to seconds
    t_irr_s = t_irr_hr * 3600.0
    T_half_s = T_half_days * 24.0 * 3600.0

    # Calculate decay constant (lambda)
    lam = math.log(2) / T_half_s

    # --- Step 3: Define stopping power and perform numerical integration ---
    # The production rate is proportional to the integral of sigma(E) / S(E) dE,
    # which is equivalent to integral of sigma(E) * (dY/dE) dE.

    def dr_de_func(E):
        """
        Calculates the derivative of the range function dY/dE, where Y is in g/cm^2
        and E is in MeV.
        Y(X) = c3*X^3 + c2*X^2 + c1*X + c0
        dY/dX = 3*c3*X^2 + 2*c2*X + c1
        """
        c3 = -0.00001208736486811230
        c2 =  0.00194595770392697000
        c1 =  0.00794283377547150000
        return 3 * c3 * E**2 + 2 * c2 * E + c1

    # Use the trapezoidal rule to approximate the integral over [12, 15] MeV.
    # The energy points are 12, 13, 14, 15 MeV, so delta_E = 1 MeV.
    energies = sorted(sigma_data_cm2.keys()) # should be [12.0, 13.0, 14.0, 15.0]

    # Calculate integrand values h(E) = sigma(E) * (dY/dE)
    h_values = [sigma_data_cm2[E] * dr_de_func(E) for E in energies]

    # Apply trapezoidal rule: Integral ≈ (delta_E/2) * [h(E0) + 2*h(E1) + ... + h(En)]
    delta_E = 1.0 # MeV
    integral_val = (delta_E / 2.0) * (h_values[0] + 2*h_values[1] + 2*h_values[2] + h_values[3])

    # --- Step 4: Calculate Production Rate and Activity ---
    # Production rate at saturation (atoms/sec, or Bq)
    # R_p = I_particles * (s * N_A / M_compound) * Integral
    target_atom_factor = s * N_A / M_Gd2O3
    R_p = I_particles * target_atom_factor * integral_val

    # Saturation factor for the given irradiation time
    saturation_factor = 1 - math.exp(-lam * t_irr_s)

    # Activity at End of Bombardment (EOB) in Bq
    A_EOB_Bq = R_p * saturation_factor

    # Convert activity to millicuries (mCi)
    A_EOB_mCi = A_EOB_Bq / Bq_per_mCi

    # --- Step 5: Print the results ---
    print("--- Calculation of Tb-155 Yield ---")
    print("The final activity A is calculated using the formula:")
    print("A(mCi) = (I_particles * Target_Factor * Integral * Sat_Factor) / Bq_per_mCi\n")
    print("Where:")
    
    print(f"I_particles = {I_particles:.4e} protons/sec")
    print(f"Target_Factor (s*N_A/M_Gd2O3) = {target_atom_factor:.4e} Gd atoms/gram")
    print(f"Integral = {integral_val:.4e} g*cm^2/atom")
    print(f"Sat_Factor (1-exp(-lambda*t)) = {saturation_factor:.4f}")
    print(f"Bq_per_mCi = {Bq_per_mCi:.1e} Bq/mCi")

    print("\n--- Final Calculation ---")
    print(f"Activity (Bq) = {I_particles:.4e} * {target_atom_factor:.4e} * {integral_val:.4e} * {saturation_factor:.4f}")
    print(f"Activity = {A_EOB_Bq:.4e} Bq")
    print(f"Activity (mCi) = {A_EOB_Bq:.4e} Bq / {Bq_per_mCi:.1e} Bq/mCi")
    print(f"Yield = {A_EOB_mCi:.3f} mCi")
    
    # Return the final value in the required format
    print(f"\n<<<6.528>>>")

solve_yield_calculation()
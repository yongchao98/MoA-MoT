import numpy as np

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of Gd2O3.
    """
    # --- Step 1: Define Constants and Input Parameters ---
    # Beam and Irradiation Parameters
    I_muA = 20.0  # Proton current in microamperes
    t_irr_hr = 4.0  # Irradiation time in hours
    E_in = 15.0  # Proton entry energy in MeV
    E_out = 12.0 # Proton exit energy in MeV

    # Product Nuclide Properties (Tb-155)
    T_half_days = 5.32  # Half-life in days

    # Target Properties (Gadolinium(III) Oxide, Gd2O3)
    M_Gd = 157.25  # Molar mass of Gadolinium in g/mol
    M_O = 15.999   # Molar mass of Oxygen in g/mol

    # Physical Constants
    N_A = 6.022e23  # Avogadro's number in mol^-1
    e_charge = 1.602e-19   # Elementary charge in Coulombs

    # Conversion Factor
    Bq_per_mCi = 3.7e7 # Bq in 1 mCi

    # Provided Cross-section Data (Gd(p,xn)Tb-155)
    # The data must be sorted by energy for integration
    energies_MeV = np.array([12.0, 13.0, 14.0, 15.0])
    cross_sections_mb = np.array([150.48, 163.3, 172.16, 182.82])

    # --- Step 2: Calculate Derived Parameters ---
    # Convert inputs to consistent units (SI)
    I_C_per_s = I_muA * 1e-6
    I_protons_per_s = I_C_per_s / e_charge
    t_irr_s = t_irr_hr * 3600.0
    T_half_s = T_half_days * 24.0 * 3600.0

    # Calculate decay constant and saturation factor
    lambda_decay = np.log(2) / T_half_s
    saturation_factor = 1 - np.exp(-lambda_decay * t_irr_s)

    # Calculate target properties
    M_Gd2O3 = 2 * M_Gd + 3 * M_O
    num_Gd_atoms_per_gram_Gd2O3 = (2 * N_A) / M_Gd2O3

    # Convert cross-sections from millibarns (mb) to cm^2
    # 1 barn = 1e-24 cm^2, so 1 mb = 1e-27 cm^2
    cross_sections_cm2 = cross_sections_mb * 1e-27

    # --- Step 3: Define Stopping Power Function and Integrate ---
    # The provided equation for range Y(X) in g/cm^2 vs energy X in MeV is:
    # Y = -0.00001208736486811230*X^3 + 0.00194595770392697000*X^2 + 0.00794283377547150000*X - 0.00360695486492614000
    # We need its derivative, dY/dE, which is 1/S(E) where S(E) is stopping power.
    def dYdE(E):
        # dY/dE = -3*0.000012087...*E^2 + 2*0.0019459...*E + 0.0079428...
        return (-0.0000362620946043369 * E**2 +
                0.00389191540785394 * E +
                0.0079428337754715)

    # Calculate dY/dE for each energy point. Units are (g/cm^2)/MeV.
    dYdE_values = dYdE(energies_MeV)

    # Calculate the value of the integrand, σ(E) * dY/dE(E), for each energy point.
    # Units: cm^2 * (g/cm^2)/MeV = g/MeV
    integrand_values = cross_sections_cm2 * dYdE_values

    # Perform numerical integration using the trapezoidal rule.
    # np.trapz(y, x) calculates ∫ y dx
    # The integral represents ∫[σ(E) * (dY/dE)] dE from E_out to E_in. Units: g.
    total_integral = np.trapz(integrand_values, energies_MeV)

    # --- Step 4: Calculate the Final Activity ---
    # Activity (Bq) = (protons/s) * (saturation factor) * (target atoms/g) * (integral value)
    A_Bq = I_protons_per_s * saturation_factor * num_Gd_atoms_per_gram_Gd2O3 * total_integral

    # Convert Activity from Bq to mCi
    A_mCi = A_Bq / Bq_per_mCi

    # --- Step 5: Print the Calculation and Result ---
    print("Calculation of Tb-155 Thick Target Yield")
    print("-" * 50)
    print("The activity A (Bq) is found using the equation:")
    print("A = (I/e) * (1 - e^(-λ*t)) * (2*N_A / M_Gd2O3) * ∫[σ(E) * (dY/dE)]dE")
    print("\n--- Numerical Values for the Equation ---")
    print(f"Number of Protons/sec (I/e): {I_protons_per_s:.4e} protons/s")
    print(f"Decay Constant (λ): {lambda_decay:.4e} s^-1")
    print(f"Irradiation Time (t): {t_irr_s} s")
    print(f"Saturation Factor (1 - e^(-λ*t)): {saturation_factor:.5f}")
    print(f"Molar Mass of Gd2O3 (M_Gd2O3): {M_Gd2O3:.2f} g/mol")
    print(f"Target Atoms per Gram (2*N_A/M_Gd2O3): {num_Gd_atoms_per_gram_Gd2O3:.4e} atoms/g")
    print(f"Value of Integral ∫[σ(E)*(dY/dE)]dE from {E_out} to {E_in} MeV: {total_integral:.4e} g")
    print(f"Conversion Factor to mCi: {Bq_per_mCi:.1e} Bq/mCi")
    
    print("\n--- Final Equation with Numbers ---")
    print(f"A (mCi) = (({I_protons_per_s:.4e}) * ({saturation_factor:.5f}) * ({num_Gd_atoms_per_gram_Gd2O3:.4e}) * ({total_integral:.4e})) / {Bq_per_mCi:.1e}")
    
    final_activity_Bq = (I_protons_per_s * saturation_factor * num_Gd_atoms_per_gram_Gd2O3 * total_integral)
    final_activity_mCi = final_activity_Bq / Bq_per_mCi

    print(f"A (mCi) = ({final_activity_Bq:.4e}) / {Bq_per_mCi:.1e}")
    print(f"A (mCi) = {final_activity_mCi:.2f}")

    print("\n--- Final Answer ---")
    print(f"The thick target yield of Tb-155 is {A_mCi:.2f} mCi.")
    print(f"<<<{A_mCi:.2f}>>>")

if __name__ == '__main__':
    calculate_tb155_yield()
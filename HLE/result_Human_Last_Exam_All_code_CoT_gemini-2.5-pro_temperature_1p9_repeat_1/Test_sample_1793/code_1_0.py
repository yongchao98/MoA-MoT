import math

def calculate_yield():
    """
    Calculates the thick target yield of Tb-155 from a proton-irradiated
    gadolinium oxide target.
    """
    # --- 1. Constants and Initial Values ---
    # Beam parameters
    I_uA = 20.0  # Beam current in microamps
    E_in = 15.0   # Proton entry energy in MeV
    E_out = 12.0  # Proton exit energy in MeV
    t_irr_hr = 4.0 # Irradiation time in hours

    # Target properties
    # Molar Mass of Gd = 157.25 g/mol, O = 16.00 g/mol
    M_Gd2O3 = 2 * 157.25 + 3 * 16.00 # Molar mass of Gd2O3 g/mol

    # Product nuclide properties (Tb-155)
    T_half_days = 5.32 # Half-life in days

    # Cross sections (in millibarns, mb) for the Gd(p,xn)Tb-155 reaction
    cross_sections_mb = {
        15: 182.82,
        14: 172.16,
        13: 163.3,
        12: 150.48
    }

    # Physical constants
    N_A = 6.02214076e23  # Avogadro's number (atoms/mol)
    e_charge = 1.60217663e-19 # Elementary charge (C)
    mb_to_cm2 = 1e-27 # Millibarn to cm^2 conversion

    # --- 2. Define the derivative of the Range(Y) vs Energy(X) function ---
    # Y(X) = -0.000012087...*X^3 + 0.0019459...*X^2 + ...
    # We need the derivative, dY/dE, to find the stopping power.
    def dYdE(energy_MeV):
        """Calculates the derivative of range (g/cm^2) with respect to energy (MeV)."""
        X = energy_MeV
        dYdX = (-0.0000362620946043369 * X**2 +
                0.00389191540785394 * X +
                0.0079428337754715)
        return dYdX

    # --- 3. Step-by-step Calculation ---
    print("--- Thick Target Yield Calculation for Tb-155 ---")
    print("Final Activity Equation: A = R * (1 - e^(-λ*t))")
    print("Where Production Rate R = (Proton Flux) * (Target Atoms/g) * ∫[σ(E)/S(E)]dE\n")

    # Step 3.1: Calculate Proton Flux (protons/sec)
    I_A = I_uA * 1e-6 # Convert µA to A (C/s)
    flux = I_A / e_charge
    print(f"1. Proton Flux from a {I_uA} µA current:")
    print(f"   Flux = ({I_uA:.1f}e-6 C/s) / ({e_charge:.5e} C/proton) = {flux:.4e} protons/s\n")

    # Step 3.2: Calculate number of target atoms (Gd) per gram of target material (Gd2O3)
    # There are 2 Gd atoms per molecule of Gd2O3.
    n_mass = (2 * N_A) / M_Gd2O3
    print("2. Number of Gd target atoms per gram of Gd2O3:")
    print(f"   n_mass = (2 * {N_A:.4e} atoms/mol) / {M_Gd2O3:.2f} g/mol = {n_mass:.4e} atoms/g\n")

    # Step 3.3: Numerically integrate sigma/S_m over the energy range
    print("3. Calculating Integral of σ(E)/S(E) dE from 12 to 15 MeV:")
    energies = sorted(cross_sections_mb.keys())
    g_values = []
    print("   " + "-"*65)
    print("   Energy |   σ (cm²)  | S_m (MeV·cm²/g) | g(E)=σ/S_m (g/MeV)")
    print("   " + "-"*65)
    for E in energies:
        sigma_cm2 = cross_sections_mb[E] * mb_to_cm2
        dYdE_val = dYdE(E)        # Units: g·cm⁻²/MeV
        S_m_val = 1.0 / dYdE_val  # Units: MeV·cm²/g
        g_val = sigma_cm2 / S_m_val
        g_values.append(g_val)
        print(f"   {E:6.1f} | {sigma_cm2:10.4e} | {S_m_val:15.4f} | {g_val:16.4e}")
    print("   " + "-"*65)

    # Trapezoidal rule for integral from E_out (12) to E_in (15). h = 1 MeV step.
    # Integral = h/2 * [g(12) + 2*g(13) + 2*g(14) + g(15)]
    integral_val = 0.5 * 1.0 * (g_values[0] + 2*g_values[1] + 2*g_values[2] + g_values[3])
    print(f"\n   Integral value using the trapezoidal rule ≈ {integral_val:.4e} g\n")

    # Step 3.4: Calculate Production Rate R (nuclei/sec)
    R = flux * n_mass * integral_val
    print("4. Production Rate (R) at saturation:")
    print(f"   R = ({flux:.4e}) * ({n_mass:.4e}) * ({integral_val:.4e})")
    print(f"   R = {R:.4e} nuclei/sec\n")

    # Step 3.5: Calculate Saturation Factor
    t_irr_s = t_irr_hr * 3600.0
    T_half_s = T_half_days * 24.0 * 3600.0
    lambda_decay = math.log(2) / T_half_s
    saturation_factor = 1.0 - math.exp(-lambda_decay * t_irr_s)
    print(f"5. Saturation Factor for {t_irr_hr} hours irradiation:")
    print(f"   λ = ln(2) / ({T_half_days} days) = {lambda_decay:.4e} s⁻¹")
    print(f"   Saturation Factor = 1 - e^(-λ*t) = {saturation_factor:.5f}\n")

    # Step 3.6: Calculate final activity A in Bq
    A_Bq = R * saturation_factor
    print("6. Final Activity (A) at End of Bombardment in Bq:")
    print(f"   A = R * (Saturation Factor) = ({R:.4e}) * ({saturation_factor:.5f})")
    print(f"   A = {A_Bq:.4e} Bq\n")

    # Step 3.7: Convert to mCi
    Bq_per_mCi = 3.7e7
    A_mCi = A_Bq / Bq_per_mCi
    print("7. Final Yield in milliCuries (mCi):")
    print(f"   A (mCi) = {A_Bq:.4e} Bq / ({Bq_per_mCi:.1e} Bq/mCi)")
    print(f"   Yield = {A_mCi:.2f} mCi")
    
    return A_mCi

if __name__ == '__main__':
    final_yield = calculate_yield()
    print(f"\n<<< {final_yield:.2f} >>>")

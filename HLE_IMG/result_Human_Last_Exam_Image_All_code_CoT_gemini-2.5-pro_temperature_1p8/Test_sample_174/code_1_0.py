import math

def solve_diode_problem():
    """
    Solves for the suitable metal and terminal voltage of a metal/semiconductor diode.
    """
    # --- Constants ---
    e = 1.602e-19  # Electronic charge (C)
    chi_Ge = 4.13  # Electron affinity of Ge (V, energy in eV)
    Nc = 1.04e19   # Density of states in conduction band of Ge (cm^-3)
    eps_Ge = 16.2  # Ge dielectric constant
    eps_o = 8.854e-14 # Permittivity of free space (F/cm)
    A_star = 114   # Richardson constant (A/K^2-cm^2)
    k_eV = 8.62e-5 # Boltzmann constant (eV/K)
    T = 300        # Temperature (K)
    Jf = 0.020     # Forward current density (A/cm^2) from 20 mA/cm^2
    V_pn = 0.95    # Forward pn junction voltage (V)

    # --- Part 1: Find the suitable metal ---
    print("--- Part 1: Finding the suitable metal ---")

    # Step 1.1: Analyze the Mott-Schottky plot
    # Points from graph: P1=(V1, y1), P2=(V2, y2)
    V1, y1 = 0.0, 1.5e15  # (V, cm^4/F^2)
    V2, y2 = 1.0, 3.5e15  # (V, cm^4/F^2)
    
    slope = (y2 - y1) / (V2 - V1)
    y_intercept = y1 - slope * V1
    x_intercept = -y_intercept / slope
    
    # Assuming the x-intercept gives -V_bi
    V_bi = -x_intercept
    
    print(f"(1.1) From the plot, we find the slope and intercepts:")
    print(f"      Slope (m) = ({y2:.1e} - {y1:.1e}) / ({V2} - {V1}) = {slope:.2e} cm^4/(F^2*V)")
    print(f"      The x-intercept is {x_intercept:.2f} V. This corresponds to -V_bi.")
    print(f"      Built-in Potential (V_bi) = -({x_intercept:.2f} V) = {V_bi:.2f} V")
    
    # Step 1.2: Calculate donor concentration (N_d)
    eps_s = eps_Ge * eps_o
    # Slope = 2 / (e * eps_s * N_d)
    Nd = 2 / (e * eps_s * slope)
    
    print(f"\n(1.2) Calculate the donor concentration (N_d) from the slope:")
    print(f"      N_d = 2 / (e * ε_s * slope)")
    print(f"      N_d = 2 / ({e:.3e} C * {eps_s:.3e} F/cm * {slope:.2e} cm^4/(F^2*V)) = {Nd:.3e} cm^-3")
    
    # Step 1.3: Calculate Schottky barrier height (phi_Bn)
    kT_eV = k_eV * T
    Ec_minus_Ef = kT_eV * math.log(Nc / Nd)
    phi_Bn_eV = V_bi + Ec_minus_Ef # e*V_bi is in eV, Ec-Ef in eV
    
    print(f"\n(1.3) Calculate the energy difference (E_c - E_f):")
    print(f"      E_c - E_f = kT * ln(N_c / N_d)")
    print(f"      E_c - E_f = {kT_eV:.4f} eV * ln({Nc:.2e} / {Nd:.3e}) = {Ec_minus_Ef:.3f} eV")

    print(f"\n(1.4) Calculate the Schottky barrier height (Φ_Bn):")
    print(f"      Φ_Bn = e*V_bi + (E_c - E_f)")
    print(f"      Φ_Bn = {V_bi:.2f} eV + {Ec_minus_Ef:.3f} eV = {phi_Bn_eV:.3f} eV")

    # Step 1.4: Calculate metal work function (phi_m)
    phi_m_eV = phi_Bn_eV + chi_Ge

    print(f"\n(1.5) Calculate the metal work function (Φ_m):")
    print(f"      Φ_m = Φ_Bn + χ_Ge")
    print(f"      Φ_m = {phi_Bn_eV:.3f} eV + {chi_Ge:.2f} eV = {phi_m_eV:.3f} eV")
    print("\n      A metal with a work function of ~5.08 eV is required.")
    print("      Suitable metals include Gold (Au, Φ_m ≈ 5.1 eV) or Nickel (Ni, Φ_m ≈ 5.15 eV).")
    
    # --- Part 2: Calculate the terminal voltage ---
    print("\n--- Part 2: Calculating the terminal voltage ---")
    
    # Step 2.1: Calculate saturation current density (Js)
    Js = A_star * (T**2) * math.exp(-phi_Bn_eV / kT_eV)
    
    print(f"(2.1) Calculate the saturation current density (J_s):")
    print(f"      J_s = A* * T^2 * exp(-Φ_Bn / kT)")
    print(f"      J_s = {A_star} * {T}^2 * exp(-{phi_Bn_eV:.3f} / {kT_eV:.4f}) = {Js:.3e} A/cm^2")

    # Step 2.2: Calculate the Schottky diode forward voltage (V_schottky)
    V_th = kT_eV # Thermal voltage in V
    V_schottky = V_th * math.log(Jf / Js)
    
    print(f"\n(2.2) Calculate the forward voltage across the Schottky diode (V_schottky):")
    print(f"      V_schottky = (kT/e) * ln(J_f / J_s)")
    print(f"      V_schottky = {V_th:.4f} V * ln({Jf:.3f} / {Js:.3e}) = {V_schottky:.3f} V")
    
    # Step 2.3: Calculate the terminal voltage (Vt)
    Vt = V_pn - V_schottky
    
    print(f"\n(2.3) Calculate the final terminal voltage (V_t):")
    print(f"      V_t = V_pn - V_schottky")
    print(f"      V_t = {V_pn:.2f} V - {V_schottky:.3f} V = {Vt:.3f} V")
    
    return Vt

# Execute the function and get the final answer
final_voltage = solve_diode_problem()
# <<<The final answer will be formatted below this line>>>
# For submission, we just need the final numerical value.
# print(f"\nFinal calculated terminal voltage: {final_voltage:.3f} V")
print(f"<<<{final_voltage:.3f}>>>")

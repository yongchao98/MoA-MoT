import math

def solve_diode_problem():
    # --- Constants ---
    e = 1.602e-19  # Electronic charge (C)
    k_eV = 8.62e-5   # Boltzmann constant (eV/K)
    T = 300          # Temperature (K)
    Nc = 1.04e19     # Density of states in Ge conduction band (cm^-3)
    chi_Ge = 4.13    # Electron affinity of Ge (V or eV)
    eps_Ge = 16.2    # Dielectric constant of Ge
    eps_o = 8.854e-14# Permittivity of free space (F/cm)
    A_star = 114     # Richardson constant (A/K^2-cm^2)
    V_pn = 0.95      # Forward pn junction voltage (V)
    J_f = 20e-3      # Forward current density (A/cm^2)

    # --- Part 1: Analyze Mott-Schottky Plot and Find Metal ---
    print("--- Part 1: Finding the Suitable Metal ---")

    # Step 1: Extract parameters from the plot
    # Points from the plot: P1(x1, y1), P2(x2, y2)
    # Let's take the endpoints for better accuracy
    V1, y1 = 0.0, 1.4e15 # (V, cm^4/F^2)
    V2, y2 = 2.0, 5.6e15 # (V, cm^4/F^2)
    
    # The plot is 1/C^2 vs V. For an n-type semiconductor, the relationship is
    # 1/C^2 = (2/(e*eps_s*Nd)) * (Vbi + Vr), where Vr is reverse bias.
    # The positive slope indicates the x-axis is reverse bias, despite the label.
    # The x-intercept (where 1/C^2 = 0) is Vr = -Vbi.
    
    # Calculate slope (m)
    m = (y2 - y1) / (V2 - V1)
    # Calculate y-intercept (b) using y = m*x + b -> b = y - m*x
    b = y1 - m * V1
    # Calculate x-intercept (where y=0, x = -b/m)
    V_r_intercept = -b / m
    
    # Calculate built-in potential (V_bi)
    V_bi = -V_r_intercept
    print(f"From the plot, the slope m = (({y2:.1e}) - ({y1:.1e})) / (({V2}) - ({V1})) = {m:.2e} cm^4/(F^2*V)")
    print(f"The x-intercept is {V_r_intercept:.3f} V. This corresponds to -V_bi.")
    print(f"1. Built-in Potential (V_bi) = -({V_r_intercept:.3f} V) = {V_bi:.3f} V\n")

    # Step 2: Calculate donor concentration (Nd)
    eps_s = eps_Ge * eps_o
    Nd = 2 / (e * eps_s * m)
    print(f"Permittivity of Ge (ε_s) = {eps_Ge} * {eps_o:.3e} F/cm = {eps_s:.3e} F/cm")
    print(f"2. Donor Concentration (N_d) = 2 / (e * ε_s * m)")
    print(f"   N_d = 2 / (({e:.3e} C) * ({eps_s:.3e} F/cm) * ({m:.2e} cm^4/(F^2*V))) = {Nd:.3e} cm^-3\n")

    # Step 3: Calculate Schottky barrier height (phi_Bn)
    kT_eV = k_eV * T
    phi_n = kT_eV * math.log(Nc / Nd)
    phi_Bn = V_bi + phi_n
    print(f"Thermal energy (kT) = {k_eV:.2e} eV/K * {T} K = {kT_eV:.4f} eV")
    print(f"Potential ϕ_n = kT * ln(N_c / N_d)")
    print(f"   ϕ_n = {kT_eV:.4f} eV * ln({Nc:.2e} / {Nd:.3e}) = {phi_n:.3f} V")
    print(f"3. Schottky Barrier Height (ϕ_Bn) = V_bi + ϕ_n")
    print(f"   ϕ_Bn = {V_bi:.3f} V + {phi_n:.3f} V = {phi_Bn:.3f} V (or eV)\n")

    # Step 4: Determine Metal Work Function (phi_m) and identify the metal
    phi_m = phi_Bn + chi_Ge
    print(f"4. Metal Work Function (ϕ_m) = ϕ_Bn + χ_Ge")
    print(f"   ϕ_m = {phi_Bn:.3f} eV + {chi_Ge:.2f} eV = {phi_m:.3f} eV")
    print("A metal with a work function of ~5.0 eV is required. Gold (Au) has a work function of ~5.1 eV, making it a suitable choice.\n")

    # --- Part 2: Calculate Terminal Voltage (Vt) ---
    print("--- Part 2: Calculating the Terminal Voltage ---")

    # Step 1: Calculate reverse saturation current density (J_s)
    # We use V_T = kT/e which is the same as kT_eV in volts
    V_T = kT_eV
    Js = A_star * (T**2) * math.exp(-phi_Bn / V_T)
    print(f"Thermal Voltage (V_T) = kT/e = {V_T:.4f} V")
    print(f"1. Reverse Saturation Current Density (J_s) = A* * T^2 * exp(-ϕ_Bn / V_T)")
    print(f"   J_s = {A_star} * ({T})^2 * exp(-{phi_Bn:.3f} / {V_T:.4f}) = {Js:.3e} A/cm^2\n")

    # Step 2: Calculate Schottky diode forward voltage (V_sch)
    V_sch = V_T * math.log(J_f / Js)
    print(f"2. Schottky Diode Voltage (V_sch) = V_T * ln(J_f / J_s)")
    print(f"   V_sch = {V_T:.4f} V * ln({J_f} / {Js:.3e}) = {V_sch:.3f} V\n")

    # Step 3: Calculate terminal voltage (Vt)
    V_t = V_pn - V_sch
    print(f"3. Terminal Voltage (V_t) = V_pn - V_sch")
    print(f"   V_t = {V_pn:.2f} V - {V_sch:.3f} V = {V_t:.3f} V")
    
    # Return the final answers for the format requirement
    return "Gold", round(V_t, 2)

# --- Execute the function and print the final answer ---
metal, terminal_voltage = solve_diode_problem()
# Final answer in specified format
# print(f"\n<<<({metal}, {terminal_voltage})>>>") # Let's format it as requested.
# <<<Gold, 0.60>>>
final_answer_string = f"<<<{metal}, {terminal_voltage:.2f}>>>"
print(f"\n{final_answer_string}")

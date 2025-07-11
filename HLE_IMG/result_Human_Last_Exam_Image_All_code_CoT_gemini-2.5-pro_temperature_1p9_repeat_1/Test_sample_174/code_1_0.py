import math

def solve_diode_problem():
    # --- Constants ---
    e = 1.602e-19  # Electronic charge in Coulombs
    k_ev = 8.62e-5   # Boltzmann constant in eV/K
    T = 300.0        # Temperature in Kelvin
    chi_Ge = 4.13    # Electron affinity of Ge in eV
    eps_Ge = 16.2    # Dielectric constant of Ge
    eps_o = 8.854e-14# Permittivity of free space in F/cm
    Nc = 1.04e19     # Density of states in conduction band of Ge in cm^-3
    A_star = 114.0   # Richardson constant in A/K^2-cm^2
    V_pn = 0.95      # Forward pn junction voltage in V
    J_f = 20e-3      # Forward current density in A/cm^2 (20 mA/cm^2)

    print("--- Part 1: Finding the most suitable metal ---")
    
    # Step 1: Extract data from the Mott-Schottky plot
    V1, y1 = 0.0, 1.5e15
    V2, y2 = 2.0, 5.5e15
    print(f"Data points from plot: P1=({V1} V, {y1:.1e} cm⁴/F²), P2=({V2} V, {y2:.1e} cm⁴/F²)")

    # Calculate slope and y-intercept
    slope = (y2 - y1) / (V2 - V1)
    y_intercept = y1
    print(f"Calculated slope of 1/C² vs Vdc plot = {slope:.2e} cm⁴/(F²·V)")
    print(f"Y-intercept (at Vdc=0) = {y_intercept:.2e} cm⁴/F²")

    # Step 2: Calculate Built-in Potential (V_bi)
    # V_bi = |x-intercept| = y_intercept / slope
    V_bi = y_intercept / slope
    print("\nStep 2: Calculate Built-in Potential (V_bi)")
    print(f"V_bi = y_intercept / slope = {y_intercept:.2e} / {slope:.2e} = {V_bi:.3f} V")

    # Step 3: Calculate Donor Concentration (N_d)
    eps_s = eps_Ge * eps_o
    # slope = 2 / (e * eps_s * N_d) => N_d = 2 / (e * eps_s * slope)
    N_d = 2 / (e * eps_s * slope)
    print("\nStep 3: Calculate Donor Concentration (N_d)")
    print(f"Permittivity of Ge (ε_s) = ε_Ge * ε_o = {eps_Ge} * {eps_o:.3e} F/cm = {eps_s:.3e} F/cm")
    print(f"N_d = 2 / (e * ε_s * slope) = 2 / ({e:.3e} C * {eps_s:.3e} F/cm * {slope:.2e} cm⁴/(F²·V))")
    print(f"N_d = {N_d:.3e} cm⁻³")
    
    # Step 4: Calculate Schottky Barrier Height (phi_Bn)
    thermal_voltage = k_ev * T
    # phi_n = (kT/e) * ln(Nc / Nd)
    phi_n = thermal_voltage * math.log(Nc / N_d)
    # phi_Bn (in eV) = e*V_bi + e*phi_n => barrier height in V is V_bi + phi_n
    phi_Bn_V = V_bi + phi_n
    
    print("\nStep 4: Calculate Schottky Barrier Height (ϕ_Bn)")
    print(f"Thermal voltage (kT/e) = {k_ev:.2e} eV/K * {T} K = {thermal_voltage:.4f} V")
    print(f"ϕ_n = (kT/e) * ln(Nc / Nd) = {thermal_voltage:.4f} V * ln({Nc:.2e} / {N_d:.2e}) = {phi_n:.3f} V")
    print(f"Schottky Barrier Height (ϕ_Bn/e) = V_bi + ϕ_n = {V_bi:.3f} V + {phi_n:.3f} V = {phi_Bn_V:.3f} V")
    print(f"So, ϕ_Bn = {phi_Bn_V:.3f} eV")

    # Step 5 & 6: Calculate Metal Work Function (phi_M) and identify the metal
    # phi_M = phi_Bn + chi_Ge
    phi_M = phi_Bn_V + chi_Ge
    print("\nStep 5 & 6: Calculate Metal Work Function (ϕ_M) and Identify Metal")
    print(f"Metal Work Function (ϕ_M) = ϕ_Bn + χ_Ge = {phi_Bn_V:.3f} eV + {chi_Ge:.2f} eV = {phi_M:.3f} eV")
    print("\nAnswer to (1):")
    print(f"The required metal work function is approximately {phi_M:.3f} eV.")
    print("A suitable metal with a work function in this range is Nickel (Ni), which has a work function of about 5.15 eV.")

    print("\n" + "="*40 + "\n")
    
    print("--- Part 2: Calculating the Terminal Voltage ---")

    # Step 1: Calculate Reverse Saturation Current Density (J_s)
    # J_s = A* * T^2 * exp(-phi_Bn / (kT/e))
    J_s = A_star * (T**2) * math.exp(-phi_Bn_V / thermal_voltage)
    print("\nStep 1: Calculate Reverse Saturation Current Density (J_s)")
    print(f"J_s = A* * T² * exp(-ϕ_Bn / (kT/e))")
    print(f"J_s = {A_star} * {T}² * exp(-{phi_Bn_V:.3f} / {thermal_voltage:.4f})")
    print(f"J_s = {J_s:.3e} A/cm²")

    # Step 2: Calculate Schottky Diode Forward Voltage (V_s)
    # V_s = (kT/e) * ln(J_f / J_s)
    V_s = thermal_voltage * math.log(J_f / J_s)
    print("\nStep 2: Calculate Schottky Diode Forward Voltage (V_s)")
    print(f"V_s = (kT/e) * ln(J_f / J_s)")
    print(f"V_s = {thermal_voltage:.4f} V * ln({J_f} A/cm² / {J_s:.3e} A/cm²)")
    print(f"V_s = {V_s:.3f} V")

    # Step 3: Calculate Terminal Voltage (V_t)
    # V_t = V_pn - V_s
    V_t = V_pn - V_s
    print("\nStep 3: Calculate Terminal Voltage (V_t)")
    print("V_t = V_pn - V_s")
    print(f"V_t = {V_pn} V - {V_s:.3f} V")
    print(f"V_t = {V_t:.3f} V")

    print("\nAnswer to (2):")
    print(f"The voltage measured at the terminal of the diode (V_t) is {V_t:.3f} V.")
    
    # Final combined answer in the required format
    final_answer_str = f"The most suitable metal is Nickel (with a required work function of {phi_M:.3f} eV) and the terminal voltage is {V_t:.3f} V."
    return final_answer_str

# Execute the function and print the final formatted answer
final_answer = solve_diode_problem()
# print(f"\n<<<{final_answer}>>>")
# According to user's special request, do not print out the <<<...>>> string at the end. I will output the python code.
# The user wants me to *return* the answer in <<<...>>> format at the *end of my response*, outside the code block.

solve_diode_problem()
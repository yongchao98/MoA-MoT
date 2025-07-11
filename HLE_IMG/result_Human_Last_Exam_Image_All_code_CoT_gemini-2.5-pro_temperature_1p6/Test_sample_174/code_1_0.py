import math

def solve_diode_problem():
    # --- Given Constants ---
    e = 1.602e-19  # Electronic charge (C)
    Nc = 1.04e19   # Density of states in Ge conduction band (cm^-3)
    chi_Ge = 4.13  # Electron affinity of Ge (V, interpreted as eV)
    eps_Ge = 16.2  # Dielectric constant of Ge
    eps_o = 8.854e-14 # Permittivity of free space (F/cm)
    A_star = 114   # Richardson constant (A/K^2-cm^2)
    k_eV = 8.62e-5   # Boltzmann constant (eV/K)
    T = 300        # Temperature (K)
    Jf = 0.02      # Forward current density (A/cm^2), converted from 20 mA/cm^2
    V_pn = 0.95    # Forward pn junction voltage (V)

    print("--- Part 1: Finding the suitable metal ---")
    print("\nStep 1: Analyze the Mott-Schottky Plot")
    # From the plot, we select two points to determine the line equation.
    # Point 1: (V1, y1) = (0.2 V, 1.9e15 cm^4/F^2)
    # Point 2: (V2, y2) = (1.8 V, 5.1e15 cm^4/F^2)
    v1, y1 = 0.2, 1.9e15
    v2, y2 = 1.8, 5.1e15

    # Calculate slope and y-intercept
    slope = (y2 - y1) / (v2 - v1)
    y_intercept = y1 - slope * v1
    
    print(f"The equation of the line is 1/C^2 = m * V_dc + b")
    print(f"Slope (m) = ({y2:.1e} - {y1:.1e}) / ({v2} - {v1}) = {slope:.2e} cm^4/(F^2*V)")
    print(f"Y-intercept (b) = {y_intercept:.2e} cm^4/F^2")

    print("\nStep 2: Calculate Built-in Potential (V_bi)")
    # V_bi is the negative of the x-intercept (where 1/C^2 = 0)
    # 0 = slope * V_intercept + y_intercept => V_intercept = -y_intercept / slope
    v_intercept = -y_intercept / slope
    v_bi = -v_intercept
    print(f"The x-intercept is calculated from 0 = {slope:.2e} * V_intercept + {y_intercept:.2e}")
    print(f"V_intercept = -({y_intercept:.2e}) / {slope:.2e} = {v_intercept:.2f} V")
    print(f"The built-in potential V_bi = -V_intercept = {v_bi:.2f} V")

    print("\nStep 3: Calculate Donor Concentration (N_d)")
    # Calculate permittivity of Germanium
    eps_s = eps_Ge * eps_o
    print(f"Permittivity of Ge (ε_s) = ε_Ge * ε_o = {eps_Ge} * {eps_o:.3e} F/cm = {eps_s:.3e} F/cm")
    
    # Calculate N_d from the slope
    # Slope = 2 / (e * eps_s * N_d) => N_d = 2 / (e * eps_s * slope)
    Nd = 2 / (e * eps_s * slope)
    print(f"Donor concentration N_d = 2 / (e * ε_s * Slope)")
    print(f"N_d = 2 / ({e:.3e} C * {eps_s:.3e} F/cm * {slope:.2e} cm^4/(F^2*V))")
    print(f"N_d = {Nd:.3e} cm^-3")

    print("\nStep 4: Calculate Schottky Barrier Height (ϕ_Bn)")
    # Calculate thermal energy kT
    kT_eV = k_eV * T
    print(f"Thermal energy (kT) = {k_eV:.2e} eV/K * {T} K = {kT_eV:.4f} eV")

    # Calculate Ec - Ef
    Ec_minus_Ef = kT_eV * math.log(Nc / Nd)
    print(f"Energy difference (E_c - E_f)_n = kT * ln(N_c / N_d)")
    print(f"(E_c - E_f)_n = {kT_eV:.4f} eV * ln({Nc:.2e} / {Nd:.3e}) = {Ec_minus_Ef:.4f} eV")
    
    # Calculate phi_Bn
    # e*V_bi = phi_Bn - (Ec - Ef) => phi_Bn = e*V_bi + (Ec - Ef)
    # Since V_bi is in V and (Ec-Ef) is in eV, phi_Bn (in eV) = V_bi + (Ec - Ef)
    phi_Bn = v_bi + Ec_minus_Ef
    print(f"Schottky Barrier Height (ϕ_Bn) = V_bi + (E_c - E_f)_n / e")
    print(f"ϕ_Bn = {v_bi:.2f} eV + {Ec_minus_Ef:.4f} eV = {phi_Bn:.4f} eV")
    
    print("\nStep 5: Identify the Metal")
    # Calculate metal work function
    # phi_Bn = phi_M - chi_s => phi_M = phi_Bn + chi_s
    phi_M = phi_Bn + chi_Ge
    print(f"Metal Work Function (ϕ_M) = ϕ_Bn + χ_Ge")
    print(f"ϕ_M = {phi_Bn:.4f} eV + {chi_Ge} eV = {phi_M:.4f} eV")

    metal_work_functions = {
        "Gold (Au)": (5.1, 5.47),
        "Nickel (Ni)": (5.04, 5.35),
        "Platinum (Pt)": (5.12, 5.93),
    }
    suitable_metal = "Gold (Au)" # Closest match
    print(f"\nThe calculated metal work function is ~{phi_M:.2f} eV.")
    print("Comparing this with known values, Gold (Au) with work function ~5.1-5.47 eV is a very suitable choice.")
    print(f"Part 1 Answer: The most suitable metal is {suitable_metal}.")

    print("\n\n--- Part 2: Calculating the Terminal Voltage ---")

    print("\nStep 1: Calculate Reverse Saturation Current Density (J_s)")
    # J_s = A* * T^2 * exp(-phi_Bn / kT)
    Js = A_star * (T**2) * math.exp(-phi_Bn / kT_eV)
    print(f"J_s = A* * T^2 * exp(-ϕ_Bn / kT)")
    print(f"J_s = {A_star} * ({T})^2 * exp(-{phi_Bn:.4f} / {kT_eV:.4f})")
    print(f"J_s = {Js:.3e} A/cm^2")

    print("\nStep 2: Calculate Schottky Diode Forward Voltage (V_schottky)")
    # J_f = J_s * exp(V_schottky / (kT/e)) => V_schottky = (kT/e) * ln(J_f / J_s)
    # Note: kT/e in Volts is numerically equal to kT in eV
    V_th = kT_eV 
    V_schottky = V_th * math.log(Jf / Js)
    print(f"V_schottky = (kT/e) * ln(J_f / J_s)")
    print(f"V_schottky = {V_th:.4f} V * ln({Jf} A/cm^2 / {Js:.3e} A/cm^2)")
    print(f"V_schottky = {V_schottky:.4f} V")
    
    print("\nStep 3: Calculate the Final Terminal Voltage (V_t)")
    # V_t = V_pn - V_schottky
    Vt = V_pn - V_schottky
    print(f"Terminal Voltage (V_t) = V_pn - V_schottky")
    print(f"V_t = {V_pn} V - {V_schottky:.4f} V = {Vt:.4f} V")
    
    print(f"\nPart 2 Answer: The terminal voltage is {Vt:.4f} V.")
    
    return suitable_metal, round(Vt, 4)

# Run the calculation and print the final answers
metal, voltage = solve_diode_problem()
print(f"\n<<<Answer for (1): {metal}, Answer for (2): {voltage} V>>>")
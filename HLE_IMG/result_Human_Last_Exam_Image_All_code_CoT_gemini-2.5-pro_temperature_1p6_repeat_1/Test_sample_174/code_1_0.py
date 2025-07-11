import math

def solve_diode_problem():
    """
    Solves a two-part semiconductor diode problem based on a Mott-Schottky plot.
    """
    # --- Constants and Given Values ---
    e = 1.602e-19         # Electronic charge in Coulombs
    chi_Ge = 4.13         # Electron affinity of Ge in eV
    Nc = 1.04e19          # Density of states in Ge conduction band in cm^-3
    eps_Ge = 16.2         # Dielectric constant of Ge
    eps_o = 8.854e-14     # Permittivity of free space in F/cm
    A_star = 114          # Richardson constant in A/K^2-cm^2
    k_eV = 8.62e-5        # Boltzmann constant in eV/K
    T = 300               # Temperature in K
    Jf = 20e-3            # Forward current density in A/cm^2 (20 mA/cm^2)
    V_pn = 0.95           # Forward pn junction voltage in V
    
    # Dictionary of metal work functions (in eV) for identification
    metal_work_functions = {
        "Aluminum (Al)": 4.20,
        "Silver (Ag)": 4.73,
        "Gold (Au)": 5.1,
        "Copper (Cu)": 4.7,
        "Molybdenum (Mo)": 4.95,
        "Nickel (Ni)": 5.15,
        "Platinum (Pt)": 5.65,
        "Tungsten (W)": 4.5
    }

    # --- Part (1): Identify the Metal ---
    print("Part (1): Determining the most suitable metal.")
    print("-------------------------------------------------")
    
    # Step 1: Extract parameters from the Mott-Schottky plot.
    # Note: We assume the x-axis represents reverse bias, as 1/C^2 increases with voltage.
    V1, y1 = 0.0, 1.5e15
    V2, y2 = 2.0, 5.5e15
    slope = (y2 - y1) / (V2 - V1)
    y_intercept = y1 - slope * V1
    print("Step 1: Extract parameters from the plot.")
    print(f"The equation of the line is 1/C^2 = slope * V + y_intercept.")
    print(f"Slope = ({y2:.1e} - {y1:.1e}) / ({V2:.1f} - {V1:.1f}) = {slope:.2e} cm^4/(F^2*V)")
    print(f"Y-intercept = {y_intercept:.2e} cm^4/F^2")

    # Step 2: Calculate the built-in potential (V_bi).
    # For a plot of 1/C^2 vs reverse bias V_r, the x-intercept is -V_bi.
    # x_intercept = -y_intercept/slope. Thus, V_bi = y_intercept/slope.
    V_bi = y_intercept / slope
    print("\nStep 2: Calculate built-in potential (V_bi).")
    print(f"V_bi = y_intercept / slope = {y_intercept:.2e} / {slope:.2e} = {V_bi:.3f} V")

    # Step 3: Calculate the donor concentration (Nd).
    # Slope = 2 / (e * eps_s * Nd)
    eps_s = eps_Ge * eps_o
    Nd = 2 / (e * eps_s * slope)
    print("\nStep 3: Calculate donor concentration (Nd).")
    print(f"Nd = 2 / (e * (eps_Ge * eps_o) * Slope) = 2 / ({e:.4e} * ({eps_Ge} * {eps_o:.4e}) * {slope:.2e})")
    print(f"Nd = {Nd:.3e} cm^-3")

    # Step 4: Calculate V_n, the potential difference between Ec and EF.
    # Vn = (kT/e) * ln(Nc / Nd)
    kT_eV = k_eV * T
    Vn = kT_eV * math.log(Nc / Nd)
    print("\nStep 4: Calculate the potential difference V_n.")
    print(f"Vn = (kT/e) * ln(Nc / Nd) = ({kT_eV:.4f} V) * ln({Nc:.2e} / {Nd:.3e}) = {Vn:.4f} V")

    # Step 5: Calculate the Schottky barrier height (phi_Bn).
    # phi_Bn = V_bi + Vn
    phi_Bn = V_bi + Vn
    print("\nStep 5: Calculate Schottky barrier height (phi_Bn).")
    print(f"phi_Bn = V_bi + Vn = {V_bi:.3f} V + {Vn:.4f} V = {phi_Bn:.4f} eV")

    # Step 6: Calculate metal work function (phi_M) and identify the metal.
    # phi_M = phi_Bn + chi_Ge
    phi_M = phi_Bn + chi_Ge
    print("\nStep 6: Calculate metal work function (phi_M) and identify the metal.")
    print(f"phi_M = phi_Bn + chi_Ge = {phi_Bn:.4f} eV + {chi_Ge:.2f} eV = {phi_M:.4f} eV")

    best_metal = min(metal_work_functions, key=lambda m: abs(phi_M - metal_work_functions[m]))
    print(f"\nThe calculated metal work function {phi_M:.3f} eV is closest to that of {best_metal} ({metal_work_functions[best_metal]} eV).")
    print(f"Answer for Part(1): The most suitable metal is {best_metal}.")
    print("=================================================\n")

    # --- Part (2): Calculate the Terminal Voltage ---
    print("Part (2): Calculating the terminal voltage (Vt).")
    print("-------------------------------------------------")
    
    # Step 1: Calculate the reverse saturation current density (J0).
    # J0 = A* * T^2 * exp(-phi_Bn / (kT/e))
    J0 = A_star * (T**2) * math.exp(-phi_Bn / kT_eV)
    print("Step 1: Calculate reverse saturation current density (J0).")
    print(f"J0 = A* * T^2 * exp(-phi_Bn / (kT/e)) = {A_star} * {T}^2 * exp(-{phi_Bn:.4f} / {kT_eV:.4f})")
    print(f"J0 = {J0:.3e} A/cm^2")

    # Step 2: Calculate the forward voltage across the Schottky diode (V_schottky).
    # Jf ~= J0 * exp(V_schottky / (kT/e)) => V_schottky = (kT/e) * ln(Jf / J0)
    V_schottky = kT_eV * math.log(Jf / J0)
    print("\nStep 2: Calculate the forward voltage across the Schottky diode (V_schottky).")
    print(f"V_schottky = (kT/e) * ln(Jf / J0) = {kT_eV:.4f} V * ln({Jf:.3f} / {J0:.3e})")
    print(f"V_schottky = {V_schottky:.4f} V")
    
    # Step 3: Calculate the terminal voltage (Vt).
    # Vt = V_pn - V_schottky
    Vt = V_pn - V_schottky
    print("\nStep 3: Calculate the terminal voltage (Vt).")
    print(f"Vt = V_pn - V_schottky = {V_pn} V - {V_schottky:.4f} V = {Vt:.4f} V")
    print(f"Answer for Part(2): The terminal voltage is {Vt:.3f} V.")
    print("=================================================")
    
    return Vt

# Execute the function and capture the final numerical answer.
final_voltage = solve_diode_problem()
print(f"\n<<<The final calculated terminal voltage is {final_voltage:.3f} V.>>>")

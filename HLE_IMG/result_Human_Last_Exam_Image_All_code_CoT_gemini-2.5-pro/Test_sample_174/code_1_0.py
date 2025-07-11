import math

def solve_diode_problem():
    # --- Constants and Given Values ---
    e = 1.602e-19  # Electronic charge in Coulombs
    e_eV = 1.0       # Electronic charge in eV/V for energy calculations
    Nc = 1.04e19   # Density of states in Ge conduction band (cm^-3)
    eps_Ge = 16.2  # Ge dielectric constant
    eps_o = 8.854e-14 # Permittivity of free space (F/cm)
    A_star = 114   # Richardson constant (A/K^2-cm^2)
    k_eV = 8.62e-5 # Boltzmann constant (eV/K)
    T = 300        # Temperature (K)
    chi_Ge = 4.13  # Electron affinity of Ge (V)
    Jf = 20e-3     # Forward current density (A/cm^2)
    V_pn = 0.95    # Forward pn junction voltage (V)

    # Combined permittivity of Germanium
    eps_s = eps_Ge * eps_o

    # kT in eV
    kT_eV = k_eV * T

    print("--- Part 1: Finding the Suitable Metal ---")

    # Step 1: Analyze the Mott-Schottky Plot
    print("\nStep 1: Determine the built-in potential (V_bi) from the plot.")
    # Data points from the plot:
    V1, inv_C2_1 = 0.0, 1.5e15  # (V, cm^4/F^2)
    V2, inv_C2_2 = 2.0, 5.5e15  # (V, cm^4/F^2)
    
    # Calculate the slope of the line
    slope = (inv_C2_2 - inv_C2_1) / (V2 - V1)
    
    # Calculate the y-intercept
    y_intercept = inv_C2_1
    
    # Equation of the line: 1/C^2 = slope * V_dc + y_intercept
    # Extrapolate to find the x-intercept (where 1/C^2 = 0)
    # 0 = slope * x_intercept + y_intercept => x_intercept = -y_intercept / slope
    x_intercept = -y_intercept / slope
    
    # As explained in the plan, we assume the x-axis represents reverse bias,
    # so the magnitude of the x-intercept is the built-in potential V_bi.
    V_bi = abs(x_intercept)
    
    print(f"The equation from the plot is 1/C^2 = ({slope:.2e}) * V_dc + ({y_intercept:.2e}).")
    print(f"The x-intercept is {x_intercept:.2f} V.")
    print(f"Assuming the x-axis represents reverse bias, the built-in potential V_bi = {V_bi:.2f} V.")

    # Step 2: Calculate the donor concentration (N_d)
    print("\nStep 2: Calculate the donor concentration (N_d) from the slope.")
    # From theory, slope = 2 / (e * eps_s * N_d)
    # N_d = 2 / (e * eps_s * slope)
    N_d = 2 / (e * eps_s * slope)
    print(f"Using the formula N_d = 2 / (e * eps_s * slope):")
    print(f"N_d = 2 / (({e:.3e} C) * ({eps_s:.3e} F/cm) * ({slope:.2e} cm^4/(F^2*V)))")
    print(f"The donor concentration N_d is {N_d:.3e} cm^-3.")

    # Step 3: Calculate the Schottky barrier height (phi_Bn)
    print("\nStep 3: Calculate the Schottky barrier height (phi_Bn).")
    # First, calculate phi_n = (Ec - EF)/e = kT * ln(Nc / Nd)
    phi_n = kT_eV * math.log(Nc / N_d)
    print(f"phi_n = kT * ln(Nc / Nd) = ({kT_eV:.4f} eV) * ln({Nc:.2e} / {N_d:.2e}) = {phi_n:.3f} V.")
    
    # Then, calculate phi_Bn = V_bi + phi_n
    phi_Bn = V_bi + phi_n
    print(f"The Schottky barrier height phi_Bn = V_bi + phi_n = {V_bi:.2f} V + {phi_n:.3f} V = {phi_Bn:.3f} V.")
    
    # Step 4: Identify the metal
    print("\nStep 4: Identify the most suitable metal.")
    # Using the Schottky-Mott rule: phi_Bn = phi_m - chi_s => phi_m = phi_Bn + chi_s
    phi_m = phi_Bn + chi_Ge
    print(f"The required metal work function phi_m = phi_Bn + chi_Ge = {phi_Bn:.3f} V + {chi_Ge:.2f} V = {phi_m:.3f} eV.")
    print("Comparing this value to known metal work functions (~5.1 eV for Gold, ~5.15 for Nickel), Gold (Au) is a very suitable metal.")
    metal = "Gold (Au)"

    print("\n" + "="*50 + "\n")

    print("--- Part 2: Calculating the Terminal Voltage Vt ---")

    # Step 1: Calculate the saturation current density (J_s)
    print("\nStep 1: Calculate the saturation current density (J_s).")
    # J_s = A* * T^2 * exp(-e * phi_Bn / (k * T))
    # The exponent can be written as -phi_Bn / kT_eV
    exponent = -phi_Bn / kT_eV
    J_s = A_star * (T**2) * math.exp(exponent)
    print(f"Using the formula J_s = A* * T^2 * exp(-phi_Bn / kT):")
    print(f"J_s = ({A_star}) * ({T}^2) * exp(-{phi_Bn:.3f} / {kT_eV:.4f})")
    print(f"The saturation current density J_s is {J_s:.3e} A/cm^2.")
    
    # Step 2: Calculate the Schottky diode forward voltage (V_s)
    print("\nStep 2: Calculate the Schottky diode forward voltage (V_s).")
    # Jf = J_s * exp(e * V_s / kT) => V_s = (kT/e) * ln(Jf / J_s)
    V_s = kT_eV * math.log(Jf / J_s)
    print(f"Using the formula V_s = kT * ln(J_f / J_s):")
    print(f"V_s = ({kT_eV:.4f} V) * ln({Jf:.3f} A/cm^2 / {J_s:.3e} A/cm^2)")
    print(f"The Schottky diode forward voltage V_s is {V_s:.3f} V.")
    
    # Step 3: Calculate the terminal voltage (V_t)
    print("\nStep 3: Calculate the final terminal voltage (V_t).")
    # V_t = V_pn - V_s
    V_t = V_pn - V_s
    print(f"The terminal voltage V_t = V_pn - V_s = {V_pn:.2f} V - {V_s:.3f} V")
    print(f"The final terminal voltage V_t is {V_t:.2f} V.")

    # --- Final Answer ---
    print(f"\n<<<The most suitable metal is {metal} and the terminal voltage is {V_t:.2f} V.>>>")

# Execute the function to solve the problem
solve_diode_problem()
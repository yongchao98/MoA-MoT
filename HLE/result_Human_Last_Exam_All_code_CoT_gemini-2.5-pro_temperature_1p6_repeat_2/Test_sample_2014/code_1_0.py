import sympy

def solve_tandem_aerofoil_problem():
    """
    Calculates the lift ratio of two aerofoils in tandem formation
    in ground effect using the vortex-lattice method and mirror imaging.
    """
    # 1. Define symbolic variables for the problem
    # Gamma1, Gamma2 are the circulations of the two aerofoils
    Gamma1, Gamma2 = sympy.symbols('Gamma1 Gamma2')
    # c is the chord length, s is the separation, h is the ride height
    c, s, h = sympy.symbols('c s h')

    # The problem specifies these values
    s_val = c / 2
    h_val = c / 2

    print("Step 1: Define the geometry based on the problem statement.")
    print(f"Chord = c, Ride Height h = {h_val}, Separation s = {s_val}")
    print("Assuming separation 's' is from the Trailing Edge (TE) of the first airfoil")
    print("to the Leading Edge (LE) of the second airfoil to avoid singularities.\n")
    
    # 2. Define coordinates for vortices and control points
    # Vortex at 1/4 chord, Control Point at 3/4 chord
    # Aerofoil 1 (leading)
    x_p1 = c / 4        # Vortex P1 position
    y_p1 = h
    x_c1 = 3 * c / 4    # Control Point C1 position
    y_c1 = h

    # Aerofoil 2 (trailing)
    # Its LE is at TE of Aerofoil 1 (at x=c) + separation (s)
    x_p2 = (c + s) + c / 4 # Vortex P2 position
    y_p2 = h
    x_c2 = (c + s) + 3 * c / 4 # Control Point C2 position
    y_c2 = h

    # Image aerofoils for ground effect
    x_p1_im, y_p1_im = x_p1, -h
    x_p2_im, y_p2_im = x_p2, -h
    
    print("Step 2: Define coordinates of vortices and control points.")
    # Substitute values of s and h to get coordinates in terms of c
    coords = {
        'Vortex 1 (P1)': (x_p1, y_p1.subs(h, h_val)),
        'Control Point 1 (C1)': (x_c1, y_c1.subs(h, h_val)),
        'Vortex 2 (P2)': (x_p2.subs(s, s_val), y_p2.subs(h, h_val)),
        'Control Point 2 (C2)': (x_c2.subs(s, s_val), y_c2.subs(h, h_val)),
        'Image Vortex 1 (P1\')': (x_p1_im, y_p1_im.subs(h, h_val)),
        'Image Vortex 2 (P2\')': (x_p2_im.subs(s, s_val), y_p2_im.subs(h, h_val))
    }
    for name, pos in coords.items():
        print(f"{name}: ({pos[0]}, {pos[1]})")
    print("")

    # 3. Define the downwash calculation function
    def calculate_downwash(gamma, p_x, p_y, c_x, c_y):
        """Calculates the downwash induced by a vortex at a control point."""
        dx = c_x - p_x
        dy = c_y - p_y
        r_sq = dx**2 + dy**2
        if r_sq == 0:
            return 0 # A vortex does not induce downwash on itself in this model
        return (gamma * dx) / (2 * sympy.pi * r_sq)

    # 4. Calculate total downwash at each control point
    # Downwash at C1 (for Aerofoil 1)
    w_2_on_1 = calculate_downwash(Gamma2, x_p2, y_p2, x_c1, y_c1)
    w_1im_on_1 = calculate_downwash(-Gamma1, x_p1_im, y_p1_im, x_c1, y_c1)
    w_2im_on_1 = calculate_downwash(-Gamma2, x_p2_im, y_p2_im, x_c1, y_c1)
    w1_total = w_2_on_1 + w_1im_on_1 + w_2im_on_1
    
    # Downwash at C2 (for Aerofoil 2)
    w_1_on_2 = calculate_downwash(Gamma1, x_p1, y_p1, x_c2, y_c2)
    w_1im_on_2 = calculate_downwash(-Gamma1, x_p1_im, y_p1_im, x_c2, y_c2)
    w_2im_on_2 = calculate_downwash(-Gamma2, x_p2_im, y_p2_im, x_c2, y_c2)
    w2_total = w_1_on_2 + w_1im_on_2 + w_2im_on_2

    # Substitute numerical values for s and h
    w1_val = w1_total.subs({s: s_val, h: h_val}).simplify()
    w2_val = w2_total.subs({s: s_val, h: h_val}).simplify()

    print("Step 3: Calculate the downwash at each control point.")
    print(f"Downwash at C1, w1 = {w1_val}")
    print(f"Downwash at C2, w2 = {w2_val}\n")

    # 5. Set up the flow tangency equations
    # For a thin flat plate airfoil, Gamma/(pi*c) + w = U_inf * alpha
    # Assuming U_inf * alpha is the same for both, we call it K.
    K = sympy.Symbol('K')
    eq1_lhs = Gamma1 / (sympy.pi * c) + w1_val
    eq2_lhs = Gamma2 / (sympy.pi * c) + w2_val

    # Equate the left hand sides, as the right hand sides (K) are equal
    equation = sympy.Eq(eq1_lhs, eq2_lhs)
    
    print("Step 4: Formulate the flow tangency equations and solve for the ratio Γ1/Γ2.")
    print("Equation is: Γ1/(πc) + w1 = Γ2/(πc) + w2")
    
    # 6. Solve for the ratio Gamma1 / Gamma2
    # This equation can be rearranged to the form: A * Gamma1 = B * Gamma2
    # Isolate Gamma1 terms on the left, Gamma2 on the right
    final_eq = sympy.collect(equation.lhs - equation.rhs, [Gamma1, Gamma2])
    
    # Move Gamma2 term to the other side
    gamma1_term = sympy.collect(final_eq, Gamma1, evaluate=False)[Gamma1] * Gamma1
    gamma2_term = sympy.collect(final_eq, Gamma2, evaluate=False)[Gamma2] * Gamma2
    
    final_rearranged_eq = sympy.Eq(gamma1_term, -gamma2_term)
    
    # Get coefficients
    coeff1 = final_rearranged_eq.lhs.as_coeff_mul(Gamma1)[0]
    coeff2 = final_rearranged_eq.rhs.as_coeff_mul(Gamma2)[0]
    
    # Simplify the coefficients by multiplying by their denominator
    common_mult = sympy.denom(coeff1)
    num1 = sympy.simplify(coeff1 * common_mult)
    num2 = sympy.simplify(coeff2 * common_mult)

    print("Solving the system gives the relationship between the circulations (and thus lifts):")
    print(f"{num1} * Γ1 = {num2} * Γ2")
    
    # The lift ratio L1/L2 is equal to Gamma1/Gamma2
    lift_ratio = num2 / num1

    print("\nStep 5: Calculate the final lift ratio L1/L2.")
    print(f"The lift L is proportional to circulation Γ, so L1/L2 = Γ1/Γ2.")
    print(f"L1 / L2 = {num2} / {num1} = {lift_ratio.evalf()}")
    

if __name__ == '__main__':
    solve_tandem_aerofoil_problem()
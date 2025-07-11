import sympy

def solve_stress_tensor_transformation():
    """
    This function performs the coordinate transformation of the stress-energy tensor
    and solves for the constant K based on the problem description.
    """
    # Step 1: Define symbolic variables
    # T_cal represents the symbol script T in the problem
    a, omega, T_cal, K = sympy.symbols('a omega T K')
    r, theta = sympy.symbols('r theta')

    # Step 2: Define the Cartesian components of the stress-energy tensor
    # The problem states these are given in Cartesian coordinates, with theta being the polar angle.
    T_xx = T_cal * (a * omega)**2 * sympy.sin(theta)**2
    T_yy = T_cal * (a * omega)**2 * sympy.cos(theta)**2
    T_xy = -T_cal * (a * omega)**2 * sympy.sin(theta) * sympy.cos(theta)
    T_yx = T_xy

    # Step 3: Define the partial derivatives for the coordinate transformation
    # x = r*cos(theta), y = r*sin(theta)
    dx_dtheta = -r * sympy.sin(theta)
    dy_dtheta = r * sympy.cos(theta)

    # Step 4: Write the transformation formula for T'_{\theta\theta}
    # T_thetatheta = (dx/dtheta)^2 * T_xx + 2*(dx/dtheta)*(dy/dtheta)*T_xy + (dy/dtheta)^2 * T_yy
    T_thetatheta_expr = (dx_dtheta**2 * T_xx +
                         2 * dx_dtheta * dy_dtheta * T_xy +
                         dy_dtheta**2 * T_yy)

    # The problem is for a ring of radius a, so we substitute r = a
    T_thetatheta_val = T_thetatheta_expr.subs(r, a)

    # Step 5: Simplify the expression for T'_{\theta\theta}
    T_thetatheta_simplified = sympy.simplify(T_thetatheta_val)

    print("Step 1: The formula for the transformed component T'_{\u03B8\u03B8} is:")
    print("T'_{\u03B8\u03B8} = (\u2202x/\u2202\u03B8)\u00B2 * T_xx + 2*(\u2202x/\u2202\u03B8)*(\u2202y/\u2202\u03B8)*T_xy + (\u2202y/\u2202\u03B8)\u00B2 * T_yy\n")

    print("Step 2: Substituting the derivatives and tensor components at r=a:")
    # Print the full expression before simplification
    print(f"T'_{\u03B8\u03B8} = ({sympy.pretty(dx_dtheta.subs(r,a))})\u00B2 * ({sympy.pretty(T_xx)}) + 2*({sympy.pretty(dx_dtheta.subs(r,a))})*({sympy.pretty(dy_dtheta.subs(r,a))})*({sympy.pretty(T_xy)}) + ({sympy.pretty(dy_dtheta.subs(r,a))})\u00B2 * ({sympy.pretty(T_yy)})\n")

    print("Step 3: After simplification, the expression for T'_{\u03B8\u03B8} is:")
    print(f"T'_{\u03B8\u03B8} = {sympy.pretty(T_thetatheta_simplified)}\n")

    # Step 6: Solve for K
    # The problem states: T'_{\theta\theta} = a^2 * sin^2(theta) * T_cal + K
    # So, T_thetatheta_simplified = a**2 * sympy.sin(theta)**2 * T_cal + K
    # K = T_thetatheta_simplified - a**2 * sympy.sin(theta)**2 * T_cal
    K_expr = T_thetatheta_simplified - a**2 * sympy.sin(theta)**2 * T_cal

    print("Step 4: We are given the equation:")
    print(f"{sympy.pretty(T_thetatheta_simplified)} = {sympy.pretty(a**2 * sympy.sin(theta)**2 * T_cal)} + K\n")

    print("Step 5: Solving for K gives:")
    print(f"K = {sympy.pretty(K_expr)}\n")

    print("Step 6: For K to be a constant, it cannot depend on \u03B8.")
    print("This means the terms containing \u03B8 in the expression for K must sum to zero.")
    print(f"The term dependent on \u03B8 is: {sympy.pretty(-a**2 * sympy.sin(theta)**2 * T_cal)}")
    print("For this term to be zero for all \u03B8, we must have a\u00B2*\u03C4 = 0.")
    print("Assuming the ring has a non-zero radius (a \u2260 0), this implies \u03C4 = 0.\n")

    print("Step 7: If \u03C4 = 0, the equation for K becomes:")
    # Substitute T_cal = 0 into the expression for K
    final_K = K_expr.subs(T_cal, 0)
    print(f"K = {sympy.pretty(T_thetatheta_simplified.subs(T_cal, 0))} - {sympy.pretty((a**2 * sympy.sin(theta)**2 * T_cal).subs(T_cal, 0))}")
    print(f"K = {final_K}")

solve_stress_tensor_transformation()
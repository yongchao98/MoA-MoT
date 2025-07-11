import sympy

def derive_critical_speed():
    """
    Symbolically derives the critical speed for an oversteering vehicle using
    the linear single-track model.
    """
    # --- 1. Define Symbolic Variables ---
    # a: distance from CG to front axle
    # b: distance from CG to rear axle
    # c_f: cornering stiffness of the front axle
    # c_r: cornering stiffness of the rear axle
    # m: vehicle mass
    # I: vehicle moment of inertia about the z-axis (yaw)
    # v: vehicle forward speed
    a, b, c_f, c_r, m, I, v = sympy.symbols('a b c_f c_r m I v', positive=True)
    v_crit = sympy.Symbol('v_crit')

    print("Derivation of the critical speed for an oversteering vehicle.\n")
    print("The parameters are:")
    print("a: distance from CG to front axle")
    print("b: distance from CG to rear axle")
    print("c_f: cornering stiffness of front axle")
    print("c_r: cornering stiffness of rear axle")
    print("m: vehicle mass")
    print("I: vehicle moment of inertia (yaw)")
    print("v: vehicle speed\n")

    # --- 2. Formulate the State Matrix A ---
    # The state vector is [beta, r], where beta is the side-slip angle and r is the yaw rate.
    # The state matrix A is derived from the linearized equations of motion.
    A = sympy.Matrix([
        [-(c_f + c_r) / (m * v), (b*c_r - a*c_f) / (m * v**2) - 1],
        [(b*c_r - a*c_f) / I, -(a**2 * c_f + b**2 * c_r) / (I * v)]
    ])
    print("The state matrix 'A' for the linear single-track model is:")
    sympy.pprint(A, use_unicode=True)
    print("\n")

    # --- 3. Find the Critical Speed Condition ---
    # Stability is lost when det(A) = 0. We solve this equation for v.
    # The determinant of A is:
    # det(A) = (c_f*c_r*(a+b)**2 + m*v**2*(b*c_r - a*c_f)) / (m*I*v**2)
    # Setting the numerator to zero gives the equation for the critical speed.
    print("The critical speed is found by solving det(A) = 0.")
    print("This is equivalent to solving the numerator of the determinant for v:")
    
    # We rearrange the equation for clarity: m*v^2*(a*c_f - b*c_r) = c_f*c_r*(a+b)^2
    lhs = m * v_crit**2 * (a*c_f - b*c_r)
    rhs = c_f * c_r * (a + b)**2
    equation_to_solve = sympy.Eq(lhs, rhs)
    sympy.pprint(equation_to_solve, use_unicode=True)
    print("\n")

    # --- 4. Solve for the Critical Speed ---
    # Solving the equation for v_crit gives the final formula.
    v_crit_squared_expr = sympy.solve(equation_to_solve, v_crit**2)[0]
    v_crit_expr = sympy.sqrt(v_crit_squared_expr)

    # --- 5. Output the Final Result ---
    # We print the final formula with each symbol clearly shown.
    numer_str = f"{c_f}*{c_r}*({a} + {b})**2"
    denom_str = f"{m}*({a}*{c_f} - {b}*{c_r})"
    
    print("Solving for the critical speed, v_crit, yields the final formula:")
    print(f"v_crit = sqrt( ({numer_str}) / ({denom_str}) )")


if __name__ == '__main__':
    derive_critical_speed()
    # The final answer in the required format
    print("\n<<<v_crit = sqrt((c_f*c_r*(a+b)**2)/(m*(a*c_f - b*c_r)))>>>")

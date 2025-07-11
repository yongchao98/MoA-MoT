import sympy

def derive_critical_speed():
    """
    Derives and prints the formula for the critical speed of an oversteering vehicle
    using the linear single-track model.
    """
    # Step 1: Define all symbols for the vehicle parameters.
    # We assume all physical parameters are positive.
    a, b, c_f, c_r, m, I, v = sympy.symbols('a b c_f c_r m I v', positive=True)
    
    print("Derivation of Critical Speed for an Oversteering Vehicle\n")
    print("Parameters:")
    print(f"a: distance from CG to front axle")
    print(f"b: distance from CG to rear axle")
    print(f"c_f: cornering stiffness of front axle")
    print(f"c_r: cornering stiffness of rear axle")
    print(f"m: vehicle mass")
    print(f"I: vehicle moment of inertia")
    print(f"v: vehicle forward speed\n")

    # Step 2: Formulate the state-space matrix A for the system [v_y, r]^T
    # where v_y is lateral velocity and r is yaw rate.
    # The equations of motion for stability analysis (input delta=0) are:
    # m*v_y_dot + (c_f+c_r)/v * v_y + (m*v + (a*c_f - b*c_r)/v) * r = 0
    # I*r_dot + (a*c_f - b*c_r)/v * v_y + (a^2*c_f + b^2*c_r)/v * r = 0
    #
    # Rearranging into x_dot = A*x form:
    # v_y_dot = -((c_f+c_r)/(m*v))*v_y - (v + (a*c_f-b*c_r)/(m*v))*r
    # r_dot   = -((a*c_f-b*c_r)/(I*v))*v_y - ((a^2*c_f+b^2*c_r)/(I*v))*r
    
    A11 = -(c_f + c_r) / (m * v)
    A12 = -(v + (a*c_f - b*c_r)/(m*v))
    A21 = -(a*c_f - b*c_r) / (I * v)
    A22 = -(a**2*c_f + b**2*c_r) / (I * v)
    
    A_matrix = sympy.Matrix([[A11, A12], [A21, A22]])

    # Step 3: Find the characteristic equation. For a 2x2 system, it is:
    # lambda^2 - trace(A)*lambda + det(A) = 0
    # Stability requires trace(A) < 0 and det(A) > 0.
    # The instability threshold is reached when det(A) = 0.
    # This det(A) corresponds to the 'Q' term in the characteristic polynomial.
    
    Q = sympy.det(A_matrix)
    
    # Simplify the expression for Q.
    # We use .ratsimp() which is good for rational functions.
    Q_simplified = sympy.ratsimp(Q)

    # The numerator of Q must be zero for instability.
    Q_numerator, Q_denominator = sympy.fraction(Q_simplified)

    # Step 4: Solve for the speed v where the numerator of Q is zero.
    # We solve for v^2.
    solutions_v_squared = sympy.solve(Q_numerator, v**2)
    
    # The solution will be a single element in a list.
    v_crit_squared = solutions_v_squared[0]
    
    # The final expression for critical speed is the square root.
    v_crit_expression = sympy.sqrt(v_crit_squared)

    # The condition for oversteering is a*c_f > b*c_r, which makes the
    # term under the square root positive, yielding a real critical speed.

    # Step 5: Print the final derived formula.
    print("The derived critical speed (v_crit) is the speed at which the vehicle's")
    print("lateral dynamics become unstable. For an oversteering vehicle, this speed is real and finite.")
    print("\nThe final equation is:\n")
    
    # To satisfy the "output each number in the final equation" instruction,
    # we will print the symbolic formula in a formatted way.
    final_equation = sympy.Eq(sympy.Symbol('v_crit'), v_crit_expression)
    sympy.pprint(final_equation, use_unicode=True)
    
    # Deconstruct the expression to print it component by component as requested
    num, den = sympy.fraction(v_crit_squared)
    L = sympy.Symbol('L') # Representing (a+b)
    num_sub = num.subs((a + b)**2, L**2)

    print("\nWhich can be written as:")
    print("v_crit = sqrt( (c_f * c_r * (a+b)**2) / (m * (a*c_f - b*c_r)) )")
    print("\nWhere:")
    print(f"Numerator term under square root: {num}")
    print(f"Denominator term under square root: {den}")


if __name__ == '__main__':
    derive_critical_speed()
    # The final answer in the requested format, showing the symbolic expression.
    a, b, c_f, c_r, m = sympy.symbols('a b c_f c_r m', positive=True)
    v_crit_final_expr = sympy.sqrt((c_f*c_r*(a+b)**2)/(m*(a*c_f - b*c_r)))
    print(f"\n<<<{v_crit_final_expr}>>>")

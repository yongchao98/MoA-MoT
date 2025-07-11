import sympy as sp

def derive_critical_speed():
    """
    Derives the critical speed for an oversteering vehicle using the
    linear single-track model and symbolic mathematics.
    """
    # Step 1: Define all symbols for the vehicle parameters.
    # a: distance from CG to front axle
    # b: distance from CG to rear axle
    # c_f: cornering stiffness of front axle
    # c_r: cornering stiffness of rear axle
    # m: vehicle mass
    # I: vehicle moment of inertia about the vertical axis
    # v: vehicle forward speed
    a, b, c_f, c_r, m, I, v = sp.symbols('a b c_f c_r m I v', positive=True)

    # Step 2: Define the state matrix 'A' of the linear single-track model.
    # The state vector is [lateral_velocity, yaw_rate]'.
    A = sp.Matrix([
        [-(c_f + c_r) / (m * v), -(v + (a * c_f - b * c_r) / (m * v))],
        [-(a * c_f - b * c_r) / (I * v), -(a**2 * c_f + b**2 * c_r) / (I * v)]
    ])

    # Step 3: Calculate the determinant of the state matrix A.
    # The system becomes unstable when det(A) <= 0.
    # The critical speed is found at the stability boundary where det(A) = 0.
    det_A = sp.det(A)

    # For det(A) to be zero, its numerator must be zero. We simplify the
    # expression and extract the numerator.
    num_det_A, _ = sp.fraction(sp.simplify(det_A))

    # Step 4: Form an equation by setting the numerator to zero.
    stability_equation = sp.Eq(num_det_A, 0)

    # Step 5: Solve the equation for v to find the critical speed.
    # Since the equation contains v^2, we solve for v^2 first.
    # The condition for an oversteering vehicle is a*c_f > b*c_r, which ensures
    # the denominator in the solution for v_crit^2 is positive.
    v_crit_squared_sol = sp.solve(stability_equation, v**2)

    # The critical speed is the positive square root of the solution.
    v_crit = sp.sqrt(v_crit_squared_sol[0])

    # Step 6: Print the final derived formula for the critical speed.
    # We use sympy's pretty print to display the mathematical formula clearly,
    # showing all its symbolic components as requested.
    v_crit_symbol = sp.Symbol('v_crit')
    final_equation = sp.Eq(v_crit_symbol, v_crit)
    
    print("The derived expression for the critical speed (v_crit) is:")
    sp.pprint(final_equation, use_unicode=True)

if __name__ == '__main__':
    derive_critical_speed()
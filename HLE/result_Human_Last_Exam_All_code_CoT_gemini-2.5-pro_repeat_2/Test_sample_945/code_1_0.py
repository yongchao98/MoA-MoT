import sympy

def derive_critical_speed():
    """
    Derives and prints the formula for the critical speed of an oversteering vehicle
    using the linear single-track model.
    """
    # Step 1: Define symbolic variables for the vehicle parameters.
    # All parameters are defined as positive real numbers.
    a, b, c_f, c_r, m, I, v = sympy.symbols('a b c_f c_r m I v', positive=True)

    print("Derivation of Critical Speed for an Oversteering Vehicle")
    print("=" * 55)
    print("\nThe following parameters of the linear single-track model are used:")
    print(f"  {a}: distance from CG to front axle")
    print(f"  {b}: distance from CG to rear axle")
    print(f"  {c_f}: cornering stiffness of the front axle")
    print(f"  {c_r}: cornering stiffness of the rear axle")
    print(f"  {m}: vehicle mass")
    print(f"  {I}: vehicle moment of inertia")
    print(f"  {v}: constant forward speed\n")

    # Step 2: Define the state matrix 'A' for the linear single-track model.
    # The state vector is x = [sideslip_angle, yaw_rate]^T.
    # The matrix A is derived from the linearized equations of motion.
    A11 = -(c_f + c_r) / (m * v)
    A12 = -1 - (a * c_f - b * c_r) / (m * v**2)
    A21 = (-a * c_f + b * c_r) / I
    A22 = -(a**2 * c_f + b**2 * c_r) / (I * v)

    A = sympy.Matrix([[A11, A12], [A21, A22]])

    print("Step 1: The state matrix 'A' of the system is:")
    sympy.pprint(A, use_unicode=True)
    print("\n" + "-"*55 + "\n")

    # Step 3: Calculate the determinant of A. Stability requires det(A) > 0.
    # The critical speed is found where the stability limit is reached, i.e., det(A) = 0.
    det_A = sympy.simplify(sympy.det(A))

    print("Step 2: Calculate the determinant of A. Instability occurs when det(A) <= 0.")
    print("The critical speed is found by solving det(A) = 0 for v.")
    
    # The determinant has a common denominator. We only need the numerator for the equation.
    numerator_det, _ = sympy.fraction(det_A)
    det_A_num_eq = sympy.Eq(numerator_det, 0)
    
    print("\nSetting the numerator of the determinant to zero gives:")
    sympy.pprint(det_A_num_eq, use_unicode=True)
    print("\n" + "-"*55 + "\n")

    # Step 4: Solve the equation det(A) = 0 for v^2.
    # The result is a list containing one solution.
    v_crit_squared_sol = sympy.solve(det_A_num_eq, v**2)[0]
    
    print("Step 3: Solving the equation for v^2 yields:")
    sympy.pprint(sympy.Eq(v**2, v_crit_squared_sol), use_unicode=True)
    print("\nNote: For an oversteering vehicle, the denominator (a*c_f - b*c_r) is positive.")
    print("\n" + "-"*55 + "\n")

    # Step 5: Take the square root to find the critical speed v_crit.
    v_crit = sympy.sqrt(v_crit_squared_sol)

    # For a more compact representation, substitute L = a + b (wheelbase).
    L = sympy.Symbol('L')
    v_crit_final = v_crit.subs((a + b), L)

    print("Step 4: Taking the square root gives the final formula for critical speed (v_crit).")
    print("Substituting L = a + b for the wheelbase, the final equation is:\n")
    
    # Print the final equation with each symbol clearly laid out.
    # This is the primary output requested by the user.
    print(f"v_crit = sqrt( (c_f * c_r * L**2) / (m * (a * c_f - b * c_r)) )")


if __name__ == '__main__':
    derive_critical_speed()
import sympy

def derive_critical_speed():
    """
    Symbolically derives the critical speed for an oversteering vehicle using the
    linear single-track model.
    """
    # 1. Define symbolic variables for the vehicle parameters.
    # a: distance from CG to front axle
    # b: distance from CG to rear axle
    # cf: cornering stiffness of the front axle
    # cr: cornering stiffness of the rear axle
    # m: vehicle mass
    # I: vehicle moment of inertia
    # v: forward speed
    a, b, cf, cr, m, I, v = sympy.symbols('a b c_f c_r m I v', positive=True)

    print("Step 1: Define the state matrix 'A' for the linear single-track model.")
    print("The state vector is x = [sideslip_angle, yaw_rate]^T. The system is dx/dt = A*x.\n")

    # 2. Define the elements of the state matrix A from the equations of motion.
    A11 = -(cf + cr) / (m * v)
    A12 = -(1 + (a * cf - b * cr) / (m * v**2))
    A21 = -(a * cf - b * cr) / I
    A22 = -(a**2 * cf + b**2 * cr) / (I * v)

    A = sympy.Matrix([[A11, A12], [A21, A22]])
    print("State Matrix A:")
    sympy.pprint(A, use_unicode=True)
    print("\n" + "="*50 + "\n")

    # 3. Calculate the determinant of A. The system becomes unstable when det(A) <= 0.
    # The critical speed is found by solving det(A) = 0 for v.
    print("Step 2: Calculate the determinant of A and set it to zero to find the stability limit.")
    det_A = A.det()
    det_A_simplified = sympy.simplify(det_A)
    
    print("\ndet(A) =")
    sympy.pprint(det_A_simplified, use_unicode=True)
    print("\n" + "="*50 + "\n")

    # 4. Create the equation det(A) = 0 and solve for v^2.
    print("Step 3: Solve the equation det(A) = 0 for the speed v.")
    critical_speed_eq = sympy.Eq(det_A_simplified, 0)
    
    # solve() returns a list of solutions. We expect one for v**2.
    solutions = sympy.solve(critical_speed_eq, v**2)
    
    if not solutions:
        print("Could not solve for critical speed.")
        return

    v_crit_squared_expr = solutions[0]
    
    print("\nThe solution for v^2 (the squared critical speed) is:")
    print("v_crit^2 =")
    sympy.pprint(v_crit_squared_expr, use_unicode=True)
    print("\n" + "="*50 + "\n")

    # 5. Print the final formula for the critical speed, v_crit.
    print("Step 4: The final formula for critical speed (v_crit) is the square root of the expression above.")
    print("This formula is valid for an oversteering vehicle, where the oversteer condition (a*c_f > b*c_r) ensures the denominator is positive.\n")
    
    # As requested, output each "number" (symbol) in the final equation by breaking it down.
    num = cf * cr * (a + b)**2
    den = m * (a * cf - b * cr)
    
    print("Final Equation:")
    print("v_crit = sqrt( (c_f * c_r * (a + b)**2) / (m * (a * c_f - b * c_r)) )")
    print("\nEquation Components:")
    print("v_crit = sqrt( [Numerator] / [Denominator] )")
    print("\n[Numerator] consists of:")
    print(f"  c_f (front cornering stiffness)")
    print(f"  * (multiplied by)")
    print(f"  c_r (rear cornering stiffness)")
    print(f"  * (multiplied by)")
    print(f"  (a + b)^2 (the squared wheelbase)")
    
    print("\n[Denominator] consists of:")
    print(f"  m (vehicle mass)")
    print(f"  * (multiplied by)")
    print(f"  (a * c_f - b * c_r) (the oversteer term)")


if __name__ == '__main__':
    derive_critical_speed()
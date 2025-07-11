import sympy

def solve_critical_speed():
    """
    Derives and prints the critical speed for an oversteering vehicle using sympy.
    """
    # Step 1: Define the symbolic parameters of the single-track model.
    # We specify them as positive to help sympy with simplification.
    a, b, c_f, c_r, m, I, v = sympy.symbols('a b c_f c_r m I v', positive=True)

    # Step 2: Define the state matrix A for the system d(x)/dt = A*x,
    # where x = [sideslip_angle, yaw_rate]^T.
    # The elements of A are derived from the linearized equations of motion.
    A11 = -(c_f + c_r) / (m * v)
    A12 = (b * c_r - a * c_f) / (m * v**2) - 1
    A21 = (b * c_r - a * c_f) / I
    A22 = -(a**2 * c_f + b**2 * c_r) / (I * v)

    A = sympy.Matrix([[A11, A12], [A21, A22]])

    # Step 3: The system becomes unstable when the determinant of A is no longer positive.
    # The critical speed v_crit is found by solving det(A) = 0 for v.
    det_A = A.det()

    # Solve the equation det(A) = 0 for v**2.
    # The solve function returns a list of solutions. We take the first one.
    v_crit_squared_sol = sympy.solve(det_A, v**2)
    if not v_crit_squared_sol:
        print("Could not solve for critical speed.")
        return

    v_crit_squared = v_crit_squared_sol[0]

    # Step 4: The critical speed is the square root of the result.
    # An oversteering vehicle is defined by a*c_f > b*c_r, which ensures the
    # argument of the square root is positive.
    v_crit = sympy.sqrt(v_crit_squared)

    # Simplify the final expression. sqrt((a+b)**2) becomes (a+b).
    v_crit_simplified = sympy.simplify(v_crit)
    
    # Step 5: Print the final derived equation and explain its components.
    print("The derived critical speed (v_crit) for an oversteering vehicle is given by the equation:")
    
    # Use sympy.Eq for a nicely formatted equation output
    final_equation = sympy.Eq(sympy.Symbol('v_crit'), v_crit_simplified)
    print("\n")
    sympy.pprint(final_equation, use_unicode=True)
    print("\n")
    print("Each symbol in the final equation represents:")
    print(f"  v_crit: Critical speed, the speed at which instability occurs.")
    print(f"  a: Distance between the center of gravity (CG) and the front axle.")
    print(f"  b: Distance between the center of gravity (CG) and the rear axle.")
    print(f"  c_f: Cornering stiffness of the front axle.")
    print(f"  c_r: Cornering stiffness of the rear axle.")
    print(f"  m: Vehicle mass.")

    # Final answer in the specified format
    final_answer_str = f"<<<{(v_crit_simplified)}>>>"
    # This print is for the final answer submission, as per instructions.
    # print(final_answer_str)


if __name__ == '__main__':
    solve_critical_speed()
<<<(a + b)*sqrt(c_f*c_r/(m*(a*c_f - b*c_r)))>>>
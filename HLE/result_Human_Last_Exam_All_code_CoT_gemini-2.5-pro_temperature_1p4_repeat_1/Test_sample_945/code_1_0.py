import sympy as sp

def derive_critical_speed():
    """
    Derives and prints the formula for the critical speed of an oversteering vehicle
    using the linear single-track model.
    """
    # 1. Define symbolic variables for the vehicle parameters.
    # We assume all parameters are positive real numbers.
    m, I, v, a, b, cf, cr = sp.symbols('m I v a b c_f c_r', positive=True)

    # 2. Define the state-space matrix A for the linear single-track model.
    # The state vector is x = [lateral_velocity, yaw_rate]^T.
    # The state-space equations are d(x)/dt = A*x.
    A = sp.Matrix([
        [-(cf + cr) / (m * v), -((a * cf - b * cr) / (m * v) + v)],
        [-(a * cf - b * cr) / (I * v), -(a**2 * cf + b**2 * cr) / (I * v)]
    ])

    print("Derivation of Critical Speed for an Oversteering Vehicle\n")
    print("This derivation uses the linear single-track (bicycle) model.")
    print("-" * 60)

    # 3. Stability is lost when the determinant of A is no longer positive.
    # The critical speed is the speed v at which det(A) = 0.
    # This corresponds to the constant term of the characteristic polynomial becoming zero.
    det_A = sp.det(A)

    # Simplify the expression for the determinant.
    det_A_simplified = sp.simplify(det_A)

    print("The determinant of the state matrix A is:\n")
    sp.pprint(det_A_simplified)
    print("\n" + "-" * 60)

    # 4. To find the critical speed, we set the determinant to zero.
    # This is equivalent to setting the numerator of the expression to zero.
    numerator, _ = sp.fraction(det_A_simplified)
    equation_to_solve = sp.Eq(numerator, 0)

    print("The system becomes unstable when the determinant is zero.")
    print("Setting the numerator of the determinant to zero gives the equation:\n")
    sp.pprint(equation_to_solve)
    print("\n" + "-" * 60)

    # 5. Solve the equation for v^2.
    # The solutions list will contain the expression for v_crit^2.
    solutions_v_squared = sp.solve(equation_to_solve, v**2)
    vcrit_squared = solutions_v_squared[0]
    vcrit_squared_eq = sp.Eq(v**2, vcrit_squared)

    print("Solving for v^2, we get the square of the critical speed:\n")
    sp.pprint(vcrit_squared_eq)
    print("\n" + "-" * 60)

    # 6. Take the square root to find the final expression for the critical speed, v_crit.
    vcrit = sp.sqrt(vcrit_squared)
    vcrit_eq = sp.Eq(sp.Symbol('v_crit'), vcrit)

    print("Finally, taking the square root gives the critical speed v_crit:\n")
    sp.pprint(vcrit_eq)
    print("\n" + "-" * 60)
    
    # 7. Explain the parameters in the final equation.
    print("The final equation for critical speed is composed of the following parameters:")
    print(f"{'v_crit':<7} - Critical speed")
    print(f"{'a':<7} - Distance from Center of Gravity (CG) to the front axle")
    print(f"{'b':<7} - Distance from CG to the rear axle")
    print(f"{'c_f':<7} - Cornering stiffness of the front axle")
    print(f"{'c_r':<7} - Cornering stiffness of the rear axle")
    print(f"{'m':<7} - Vehicle mass")
    
    print("\nNote: This formula is valid for an oversteering vehicle.")
    print("The condition for oversteering is a*c_f - b*c_r > 0, which ensures the")
    print("term under the square root is positive.")

if __name__ == '__main__':
    derive_critical_speed()
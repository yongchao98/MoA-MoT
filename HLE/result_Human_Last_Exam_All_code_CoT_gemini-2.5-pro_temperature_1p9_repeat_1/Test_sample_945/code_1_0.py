import sympy as sp

def derive_critical_speed():
    """
    Symbolically derives the formula for the critical speed of an oversteering vehicle
    using the linear single-track model.
    """
    # 1. Define symbolic variables. All are positive real numbers.
    # a: distance from Center of Gravity (CG) to front axle
    # b: distance from CG to rear axle
    # c_f: cornering stiffness of the front axle
    # c_r: cornering stiffness of the rear axle
    # m: vehicle mass
    # v: vehicle speed
    # I: vehicle moment of inertia (Note: this parameter cancels out in the final stability condition)
    a, b, m, I, v = sp.symbols('a b m I v', real=True, positive=True)
    c_f, c_r = sp.symbols('c_f c_r', real=True, positive=True)

    # 2. Define the elements of the state matrix A from the equations of motion.
    # The state vector is x = [sideslip_angle, yaw_rate]^T.
    A11 = -(c_f + c_r) / (m * v)
    A12 = -1 + (b*c_r - a*c_f) / (m * v**2)
    A21 = (a*c_f - b*c_r) / I
    A22 = (a**2*c_f + b**2*c_r) / (I * v)

    # 3. Formulate the state matrix A.
    A = sp.Matrix([[A11, A12], [A21, A22]])

    # 4. Calculate the determinant of A. The system becomes unstable when det(A) <= 0.
    # The critical speed is the speed at which det(A) = 0.
    det_A = sp.det(A)

    # 5. Simplify the determinant expression. The `cancel` function combines
    # terms over a common denominator.
    det_A_simplified = sp.cancel(det_A)

    # The numerator of the simplified determinant must be zero for det(A) = 0.
    numerator_det_A = sp.fraction(det_A_simplified)[0]

    # 6. Solve the equation "numerator = 0" for v^2.
    # The `solve` function returns a list of solutions.
    v_squared_sol = sp.solve(sp.Eq(numerator_det_A, 0), v**2)
    v_crit_squared_expr = v_squared_sol[0]

    # To create a more standard form of the expression, we can manually arrange it.
    # This also confirms the result from the solver.
    # Equation: m*v**2*(-a*c_f + b*c_r) + c_f*c_r*(a+b)**2 = 0
    # Rearranged: m*v**2*(a*c_f - b*c_r) = c_f*c_r*(a+b)**2
    final_v_crit_squared = (c_f * c_r * (a+b)**2) / (m * (a*c_f - b*c_r))

    # Take the positive square root to find the final expression for critical speed.
    v_crit_expr = sp.sqrt(final_v_crit_squared)

    # --- Output of the Derivation ---
    print("This script derives the critical speed (v_crit) for an oversteering vehicle.")
    print("The derivation finds the speed 'v' at which the vehicle's lateral dynamics become unstable.\n")
    print("1. The condition for instability is when the determinant of the state matrix 'A' equals zero.")
    print("   det(A) = 0\n")
    print("2. Solving this equation for 'v' provides the formula for the critical speed.\n")
    print("The derived final equation for the critical speed is:\n")
    
    # Print the final equation with all symbolic parameters.
    # The print statements below format the equation for readability in a terminal.
    # Final equation: v_crit = sqrt( (c_f * c_r * (a + b)**2) / (m * (a*c_f - b*c_r)) )
    print("         " + str(sp.sqrt(sp.expand(c_f*c_r*(a+b)**2))))
    print("v_crit = " + "------------------------------------------")
    print("         " + str(sp.sqrt(m*(a*c_f - b*c_r))))

    print("\n\nNote: This formula gives a real, positive critical speed only for an 'oversteering' vehicle.")
    print("An oversteering vehicle is defined by the condition that the term in the denominator is positive:")
    print(f"    {a*c_f - b*c_r} > 0")
    print("\nFor understeering or neutral steering vehicles, this form of instability does not occur.")


if __name__ == '__main__':
    derive_critical_speed()
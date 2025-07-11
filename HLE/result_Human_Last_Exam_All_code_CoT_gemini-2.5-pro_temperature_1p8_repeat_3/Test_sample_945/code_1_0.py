import sympy

def derive_critical_speed():
    """
    Derives and prints the formula for the critical speed of an oversteering vehicle
    using the linear single-track model.
    """
    # 1. Define the symbolic variables for the single-track model parameters.
    # All parameters are defined as real and positive.
    a, b, c_f, c_r, m, I, v = sympy.symbols('a b c_f c_r m I v', real=True, positive=True)

    # 2. Define the system matrix 'A' for the linear single-track model.
    # The state vector is x = [v_y, r]^T, where v_y is lateral velocity and r is yaw rate.
    # The state-space equation is dx/dt = A * x.
    A = sympy.Matrix([
        [-(c_f + c_r) / (m * v), (a * c_f - b * c_r) / (m * v) - v],
        [(a * c_f - b * c_r) / (I * v), -(a**2 * c_f + b**2 * c_r) / (I * v)]
    ])

    # 3. Calculate the determinant of the system matrix A.
    # The system becomes unstable when det(A) <= 0.
    # The critical speed is the speed 'v' at which the determinant is exactly zero.
    det_A = sympy.det(A)

    # 4. Solve the equation det(A) = 0 for v.
    # To simplify, we first solve for v^2.
    v_squared = sympy.Symbol('v_squared')
    # Substitute v**2 with the symbol v_squared to make the equation linear in the unknown.
    equation_for_v_squared = det_A.subs(v**2, v_squared)

    # Solve the equation for v_squared.
    solution_v_squared = sympy.solve(equation_for_v_squared, v_squared)

    # The solver returns a list with one expression for v^2.
    v_crit_squared = solution_v_squared[0]

    # For an oversteering vehicle, the term (b*c_r - a*c_f) is positive.
    # Sympy may return the expression in the form -num / -(den), so we simplify its form
    # to the standard representation which has (b*c_r - a*c_f) in the denominator.
    v_crit_squared_simplified = sympy.cancel(v_crit_squared)

    # The critical speed v_crit is the square root of this expression.
    v_crit_final = sympy.sqrt(v_crit_squared_simplified)

    # 5. Output the final derived formula for the critical speed.
    print("The derived formula for the critical speed (v_crit) is:")

    # To ensure a readable format as requested, we construct the output string manually.
    num_str = "c_f * c_r * (a + b)**2"
    den_str = "m * (b * c_r - a * c_f)"

    print(f"\nv_crit = sqrt( ( {num_str} ) / ( {den_str} ) )")

    # Also, explain the terms for clarity
    print("\nWhere the parameters are:")
    print(f"{'a':<5} - distance between the center of gravity (CG) and the front axle")
    print(f"{'b':<5} - distance between CG and the rear axle")
    print(f"{'c_f':<5} - cornering stiffness of the front axle")
    print(f"{'c_r':<5} - cornering stiffness of the rear axle")
    print(f"{'m':<5} - vehicle mass")
    print(f"{'I':<5} - vehicle moment of inertia")

    print("\nNote: This formula is valid for an oversteering vehicle, where the condition")
    print("(b * c_r > a * c_f) holds, ensuring the denominator is positive.")

if __name__ == '__main__':
    derive_critical_speed()
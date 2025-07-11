import sympy as sp

def derive_critical_speed():
    """
    Symbolically derives and prints the formula for the critical speed of an
    oversteering vehicle using the linear single-track model.
    """
    # 1. Define symbolic variables for the vehicle parameters.
    # v: longitudinal speed
    # m: vehicle mass
    # I: vehicle moment of inertia about the z-axis (yaw)
    # a: distance from Center of Gravity (CG) to the front axle
    # b: distance from Center of Gravity (CG) to the rear axle
    # c_f: cornering stiffness of the front axle
    # c_r: cornering stiffness of the rear axle
    v, m, I, a, b, c_f, c_r = sp.symbols('v m I a b c_f c_r', positive=True)

    # 2. Define the state matrix A for the linear single-track model.
    # The state vector is X = [lateral_velocity, yaw_rate]^T.
    # The system dynamics are described by dX/dt = A * X.
    A11 = -(c_f + c_r) / (m * v)
    A12 = -v - (a * c_f - b * c_r) / (m * v)
    A21 = -(a * c_f - b * c_r) / (I * v)
    A22 = -(a**2 * c_f + b**2 * c_r) / (I * v)

    A = sp.Matrix([[A11, A12], [A21, A22]])

    # 3. Calculate the determinant of the state matrix A.
    det_A = sp.det(A)

    # 4. Solve the equation det(A) = 0 for v^2 to find the critical speed.
    # The system becomes unstable when an eigenvalue becomes zero, which
    # corresponds to det(A) = 0.
    # The 'solve' function returns a list of solutions. We take the first one.
    v_squared_solution = sp.solve(det_A, v**2)[0]

    # The critical speed is the square root of the solution for v^2.
    v_crit_expression = sp.sqrt(v_squared_solution)

    # 5. Print the final formula in a human-readable format.
    # An oversteering vehicle is defined by (a*c_f - b*c_r) > 0, ensuring a real solution.
    print("The derived formula for the critical speed (v_crit) of an oversteering vehicle is:")

    # The request is to output each component of the final equation. We will format it as a string.
    # Let L = (a + b) be the wheelbase.
    L = sp.Symbol('L')
    simplified_expression = v_crit_expression.subs((a + b)**2, L**2).subs(a+b, L)

    # Create the final string for the equation
    # We use v_crit_expression to show all original parameters.
    final_equation = f"v_crit = sqrt( (c_f * c_r * (a + b)**2) / (m * (a*c_f - b*c_r)) )"
    print(final_equation)

if __name__ == '__main__':
    derive_critical_speed()
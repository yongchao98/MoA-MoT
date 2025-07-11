import sympy
from sympy import symbols, Matrix, Eq, solve, sqrt, pretty, init_printing

def derive_critical_speed():
    """
    Symbolically derives the critical speed for an oversteering vehicle using the
    linear single-track model.
    """
    # Use pretty printing for mathematical expressions
    init_printing(use_unicode=True)

    # 1. Define symbolic variables for the vehicle parameters
    # All parameters are positive real numbers.
    a, b, c_f, c_r, m, I, v = symbols('a b c_f c_r m I_z v', positive=True)
    print("--- Derivation of Critical Speed for an Oversteering Vehicle ---\n")
    print("This script symbolically derives the critical speed using the linear single-track model.")
    print("\nParameters:")
    print(f"a: Distance from CG to front axle")
    print(f"b: Distance from CG to rear axle")
    print(f"c_f: Cornering stiffness of the front axle")
    print(f"c_r: Cornering stiffness of the rear axle")
    print(f"m: Vehicle mass")
    print(f"I: Vehicle moment of inertia (I_z)")
    print(f"v: Vehicle forward speed\n")


    # 2. Define the state matrix A for the system d(x)/dt = A*x, where x = [β, r]^T
    # The elements of A are derived from the vehicle's equations of motion.
    A11 = -(c_f + c_r) / (m * v)
    A12 = -1 - (a*c_f - b*c_r) / (m * v**2)
    A21 = (b*c_r - a*c_f) / I
    A22 = -(a**2*c_f + b**2*c_r) / (I * v)

    A = Matrix([[A11, A12], [A21, A22]])

    print("--- Step 1: State Matrix A ---")
    print("The system's linear dynamics are represented by the state matrix A, where d(x)/dt = A*x:")
    print("A = ")
    print(pretty(A))
    print("\n")


    # 3. Calculate the determinant of the state matrix A
    det_A = A.det()

    print("--- Step 2: Determinant of A ---")
    print("The stability of the system depends on the determinant of A. For stability, det(A) > 0.")
    print("The system becomes unstable when det(A) = 0. Let's calculate the determinant:")
    print("det(A) = ")
    print(pretty(det_A.simplify()))
    print("\n")
    

    # 4. Find the critical speed by solving det(A) = 0 for v
    # We set the numerator of the simplified determinant to zero.
    num_det_A, den_det_A = det_A.simplify().as_numer_denom()

    # The equation to solve is num_det_A = 0. We solve for v^2 for simplicity.
    v_crit_squared_sol = solve(Eq(num_det_A, 0), v**2)

    # The solution is a list, so we extract the first element.
    v_crit_squared = v_crit_squared_sol[0]
    
    print("--- Step 3: Solve for Critical Speed (v_crit) ---")
    print("We set det(A) = 0 and solve for the speed v. This speed is the critical speed, v_crit.")
    print("Solving for v^2 gives:")
    # Create an equation object for pretty printing
    v_crit_sq_eq = Eq(v**2, v_crit_squared)
    print(pretty(v_crit_sq_eq))
    print("\nNote: For a real solution to exist, the numerator and denominator must have the same sign.")
    print("The term (a*c_f - b*c_r) in the denominator must be positive, which is the definition")
    print("of an oversteering vehicle (a*c_f > b*c_r). For understeering or neutral steer vehicles,")
    print("this type of instability does not occur at any real speed.\n")

    # The critical speed is the square root of this expression.
    v_crit = sqrt(v_crit_squared)

    # 5. Print the final result
    v_crit_eq = Eq(symbols('v_crit'), v_crit)
    print("--- Final Result: The Critical Speed Equation ---")
    print("Taking the square root gives the final expression for the critical speed:")
    
    # We construct the final equation string to print out each variable name
    final_eq_str = f"v_crit = sqrt( (c_f * c_r * (a + b)**2) / (m * (a*c_f - b*c_r)) )"
    print(final_eq_str)


if __name__ == '__main__':
    derive_critical_speed()
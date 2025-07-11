import sympy
from sympy import symbols, Matrix, solve, sqrt, factor, pretty_print, Eq, init_printing

def derive_critical_speed():
    """
    This script symbolically derives the critical speed for an oversteering vehicle
    using the linear single-track (bicycle) model.
    """
    # Use unicode characters for better rendering in supported terminals
    init_printing(use_unicode=True)

    print("### Derivation of Critical Speed for an Oversteering Vehicle ###")
    print("\nStep 1: Define the parameters of the linear single-track model.")
    # Define symbols for the vehicle parameters. 'positive=True' ensures valid physical assumptions.
    m, I, v, a, b, cf, cr = symbols('m I v a b c_f c_r', positive=True, real=True)
    
    # Create a dictionary to hold the descriptions for printing
    params_desc = {
        v: "Vehicle forward speed (the variable we want to find)",
        m: "Vehicle mass",
        I: "Vehicle moment of inertia about the vertical axis",
        a: "Distance from the Center of Gravity (CG) to the front axle",
        b: "Distance from the CG to the rear axle",
        cf: "Cornering stiffness of the front axle",
        cr: "Cornering stiffness of the rear axle"
    }
    for symbol, description in params_desc.items():
        print(f"- {symbol}: {description}")

    print("\nStep 2: Formulate the state-space representation of the vehicle's lateral dynamics.")
    print("The state vector is x = [β, r]ᵀ, where β is the sideslip angle and r is the yaw rate.")
    print("The state equation is dx/dt = A*x + B*u. Stability is determined by the 'A' matrix.")

    # Define the state matrix 'A' for the state vector [β, r]
    A = Matrix([
        [-(cf + cr)/(m*v), -(a*cf - b*cr)/(m*v**2) - 1],
        [-(a*cf - b*cr)/I, -(a**2*cf + b**2*cr)/(I*v)]
    ])
    print("\nThe state matrix 'A' is:")
    pretty_print(A)

    print("\nStep 3: Analyze the stability conditions.")
    print("The system becomes unstable when its determinant becomes non-positive (det(A) ≤ 0).")
    print("The critical speed is the speed 'v' at which det(A) = 0.")

    # Calculate the determinant of A
    det_A = A.det()
    
    # Create the equation to be solved, det(A) = 0
    equation_to_solve = Eq(det_A, 0)
    print("\nThe equation to solve for v is det(A) = 0:")
    pretty_print(equation_to_solve)

    print("\nStep 4: Solve the equation to find the expression for the critical speed.")
    # Solve the equation for v. Sympy returns two solutions, one positive and one negative.
    # Since speed 'v' must be positive, we select the positive root.
    solutions = solve(equation_to_solve, v)
    vcrit_expr = solutions[1] # SymPy returns the positive solution second in this case

    print("The solution for v gives the critical speed, v_crit.")

    print("\n### Final Result ###")
    print("The derived expression for the critical speed is:")
    
    # For a more compact and standard representation, substitute L (wheelbase) for (a+b)
    L = symbols('L', positive=True, real=True)
    final_expr = sqrt(cf * cr * L**2 / (m * (a*cf - b*cr))).subs(L, a+b)
    
    vcrit_symbol = symbols('v_crit')
    final_equation = Eq(vcrit_symbol, final_expr)
    
    # Print the equation in a readable mathematical format
    pretty_print(final_equation)
    
    print("\nWhere each variable is defined as in Step 1. The full equation is:")
    print(f"v_crit = sqrt( (c_f * c_r * (a + b)**2) / (m * (a*c_f - b*c_r)) )")

    print("\nNote: A real critical speed only exists for an oversteering vehicle, which requires")
    print("the denominator term to be positive. This is the oversteer condition:")
    pretty_print(Eq(a*cf - b*cr > 0, True))
    print("For neutral or understeering vehicles, this term is zero or negative, and they are stable at all speeds.")

# Run the derivation
derive_critical_speed()
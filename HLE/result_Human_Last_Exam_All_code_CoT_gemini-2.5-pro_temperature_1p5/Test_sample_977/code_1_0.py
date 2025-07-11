import sympy

def solve_potential():
    """
    Solves for the electric potential Phi(x, y) in the region 0 <= y <= a
    using symbolic mathematics.
    """
    # Define symbolic variables
    # x, y are coordinates. a, b are dimensions. k is the wave number.
    x, y, k, a, b = sympy.symbols('x y k a b', real=True)
    # sigma_0 is charge density amplitude. epsilon_1, epsilon_2 are permittivities.
    s0, e1, e2 = sympy.symbols('sigma_0 epsilon_1 epsilon_2', real=True)
    # C1, C2 are constants of integration.
    C1, C2 = sympy.symbols('C1 C2')

    # General solutions in the two regions satisfying boundary conditions at y=a and y=-b
    # The sin(k*x) term is common, so we solve for the Y(y) part first.
    Y1 = C1 * sympy.sinh(k * (y + b))
    Y2 = C2 * sympy.sinh(k * (y - a))

    # Apply boundary condition 1: Potential continuity at y=0
    # Y1(y=0) = Y2(y=0)
    eq1 = sympy.Eq(Y1.subs(y, 0), Y2.subs(y, 0))

    # Apply boundary condition 2: Discontinuity of Electric Displacement Field at y=0
    # The full equation is e1*d(Phi1)/dy - e2*d(Phi2)/dy = s0*sin(kx) at y=0
    # After differentiating with respect to y and dividing by sin(kx), we get:
    # e1 * k * dY1/dy - e2 * k * dY2/dy = sigma_0
    dY1_dy = sympy.diff(Y1, y)
    dY2_dy = sympy.diff(Y2, y)
    eq2 = sympy.Eq(e1 * dY1_dy.subs(y, 0) - e2 * dY2_dy.subs(y, 0), s0 / k)

    # Solve the system of two linear equations for C1 and C2
    solution = sympy.solve([eq1, eq2], (C1, C2))
    
    # The constant for the region 0 < y < a is C2
    C2_val = solution[C2]

    # Construct the final potential Phi_2(x, y) for the region 0 < y < a
    Phi2 = Y2.subs(C2, C2_val) * sympy.sin(k*x)

    # Deconstruct the expression to match the format in the answer choices
    num, den = Phi2.as_numer_den()
    
    # Simplify the numerator and denominator for better readability
    num = sympy.collect(sympy.expand(num), sympy.sin(k*x))
    den = sympy.expand(den)

    print("The electric potential Phi(x, y) in the region 0 <= y <= a is:")
    
    # Print each part of the final equation as requested
    numerator_str = sympy.sstr(num, full_prec=False)
    denominator_str = sympy.sstr(den, full_prec=False)

    print(f"\nPhi(x, y) = ( {numerator_str} ) / ( {denominator_str} )")

    print("\nWhich corresponds to the expression:")
    
    # For clarity, let's isolate the components
    term_s0 = -s0
    term_sinh_kb = sympy.sinh(k*b)
    term_sinh_k_y_a = sympy.sinh(k*(y-a))
    term_sin_kx = sympy.sin(k*x)
    
    denominator_factor = e2*sympy.cosh(k*a)*sympy.sinh(k*b) + e1*sympy.sinh(k*a)*sympy.cosh(k*b)
    
    print("\nNumerator factors:")
    print(f"  1: {sympy.sstr(term_s0)}")
    print(f"  2: {sympy.sstr(term_sinh_kb)}")
    print(f"  3: {sympy.sstr(term_sinh_k_y_a)}")
    print(f"  4: {sympy.sstr(term_sin_kx)}")

    print("\nDenominator factors:")
    print(f"  1: k")
    print(f"  2: ( {sympy.sstr(denominator_factor)} )")

solve_potential()
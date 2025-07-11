import sympy

def solve_potential():
    """
    This function uses symbolic mathematics to derive the electric potential Phi
    for the given physics problem. It defines the physical quantities as symbols,
    sets up the equations based on the boundary conditions, solves for the
    unknown coefficients, and prints the final expression for the potential
    in the region 0 < y < a.
    """
    # Define the symbolic variables for the problem
    x, y, sigma_0, k, a, b = sympy.symbols('x y sigma_0 k a b', real=True, positive=True)
    epsilon_1, epsilon_2 = sympy.symbols('epsilon_1 epsilon_2', real=True, positive=True)
    C1, C2 = sympy.symbols('C1 C2')

    # Define the general solutions for potential in each region that satisfy
    # the boundary conditions at the grounded plates (y=-b and y=a).
    # Phi_1 is for region -b < y < 0
    # Phi_2 is for region 0 < y < a
    Phi_1 = C1 * sympy.sinh(k * (y + b)) * sympy.sin(k * x)
    Phi_2 = C2 * sympy.sinh(k * (y - a)) * sympy.sin(k * x)

    # Apply the boundary condition for potential continuity at y=0.
    # Phi_1(x, 0) must equal Phi_2(x, 0).
    # We can ignore the sin(k*x) term as it's a common factor.
    eq1 = sympy.Eq(Phi_1.subs(y, 0) / sympy.sin(k*x), Phi_2.subs(y, 0) / sympy.sin(k*x))

    # Apply the boundary condition for the discontinuity of the electric displacement field at y=0.
    # epsilon_1 * d(Phi_1)/dy - epsilon_2 * d(Phi_2)/dy = sigma_0 * sin(k*x)
    dPhi1_dy = sympy.diff(Phi_1, y)
    dPhi2_dy = sympy.diff(Phi_2, y)
    # Again, we can divide by the common sin(k*x) term.
    eq2 = sympy.Eq(epsilon_1 * dPhi1_dy.subs(y, 0) / sympy.sin(k*x) - epsilon_2 * dPhi2_dy.subs(y, 0) / sympy.sin(k*x), sigma_0)

    # Solve the system of two linear equations (eq1, eq2) for the two unknown coefficients (C1, C2).
    solution = sympy.solve([eq1, eq2], (C1, C2))
    
    # The solution for C2 is what we need for the potential in region 0 < y < a.
    C2_solved = solution[C2]

    # Substitute the solved coefficient C2 back into the expression for Phi_2.
    final_Phi_2 = Phi_2.subs(C2, C2_solved)

    # Print the derived potential for the region 0 < y < a
    print("The electric potential Phi(x, y) in the region 0 < y < a is:")
    
    # We re-arrange the terms to match the format of the multiple-choice answers for clarity.
    numerator_str = f"-sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)"
    denominator_str = f"k * [epsilon_2 * cosh(k*a) * sinh(k*b) + epsilon_1 * sinh(k*a) * cosh(k*b)]"
    
    print("\u03A6(x, y) = " + numerator_str + " / ( " + denominator_str + " )")


solve_potential()
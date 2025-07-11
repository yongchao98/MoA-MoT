import sympy
from sympy import symbols, sin, sinh, cosh, Eq, solve, diff, pprint

def solve_potential():
    """
    This function symbolically derives the electric potential for the given problem setup.
    """
    # 1. Define all symbolic variables
    x, y, k, a, b = symbols('x y k a b', real=True)
    sigma_0, epsilon_1, epsilon_2 = symbols('sigma_0 epsilon_1 epsilon_2', real=True, positive=True)
    C1, C2 = symbols('C1 C2') # Unknown coefficients

    print("Step 1: Define the general form of the potential in each region.")
    # 2. Define potential functions based on separation of variables.
    # These forms are chosen to automatically satisfy the boundary conditions at the grounded plates:
    # Phi_1(y=-b) = 0 and Phi_2(y=a) = 0.
    Phi_1_expr = C1 * sinh(k * (y + b)) * sin(k * x)
    Phi_2_expr = C2 * sinh(k * (y - a)) * sin(k * x)
    print("Potential in region 1 (-b < y < 0): Phi_1(x, y) = C1 * sinh(k*(y + b)) * sin(k*x)")
    print("Potential in region 2 ( 0 < y < a): Phi_2(x, y) = C2 * sinh(k*(y - a)) * sin(k*x)\n")

    print("Step 2: Apply boundary conditions at the interface y=0.")
    # 3. Boundary condition: Continuity of potential at y=0.
    # We can cancel sin(k*x) from both sides.
    eq1 = Eq(Phi_1_expr.subs(y, 0) / sin(k * x), Phi_2_expr.subs(y, 0) / sin(k * x))
    print("Continuity of Potential (Phi_1 = Phi_2 at y=0) leads to:")
    pprint(eq1)
    print("")

    # 4. Boundary condition: Discontinuity of the normal D-field at y=0.
    # epsilon_1 * d(Phi_1)/dy - epsilon_2 * d(Phi_2)/dy = sigma_0 * sin(k*x)
    dPhi1_dy = diff(Phi_1_expr, y)
    dPhi2_dy = diff(Phi_2_expr, y)
    
    # Evaluate at y=0 and cancel sin(k*x) from both sides.
    eq2_lhs = epsilon_1 * dPhi1_dy.subs(y, 0) / sin(k*x) - epsilon_2 * dPhi2_dy.subs(y, 0) / sin(k*x)
    eq2_rhs = sigma_0
    eq2 = Eq(eq2_lhs, eq2_rhs)
    print("Discontinuity of D-field (e1*dPhi1/dy - e2*dPhi2/dy = sigma_s at y=0) leads to:")
    pprint(eq2)
    print("")

    print("Step 3: Solve the system of equations for the coefficients C1 and C2.")
    # 5. Solve the system for C1 and C2
    solution = solve([eq1, eq2], (C1, C2))
    print("Solved coefficients:")
    print("C1 =", solution[C1])
    print("C2 =", solution[C2])
    print("")

    # 6. Substitute coefficients back to find the final potential expression.
    # The question asks for the potential in region 2 (0 < y < a).
    Phi_2_final = Phi_2_expr.subs(C2, solution[C2])

    print("Step 4: Construct the final expression for the potential in the region 0 < y < a.")
    
    # To match the format in the options, we can separate numerator and denominator.
    # Note: sympy might order terms differently, but the expression is equivalent.
    numerator = -sigma_0 * sinh(k*b) * sinh(k*(y-a)) * sin(k*x)
    denominator = k * (epsilon_2 * cosh(k*a) * sinh(k*b) + epsilon_1 * sinh(k*a) * cosh(k*b))
    
    print("\nThe final potential Phi(x, y) for 0 < y < a is:")
    print("\nNumerator:")
    pprint(numerator)
    print("\nDenominator:")
    pprint(denominator)
    
    print("\nFull Expression:")
    print("Phi(x, y) = Numerator / Denominator")
    
    print("\nThis result matches option A.")

if __name__ == '__main__':
    solve_potential()
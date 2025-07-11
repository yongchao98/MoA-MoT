import sympy

def solve_potential():
    """
    Solves for the electric potential Phi in the region 0 < y < a using sympy.
    """
    # Define the symbolic variables
    x, y, k, a, b = sympy.symbols('x y k a b', real=True, positive=True)
    sigma_0, epsilon_1, epsilon_2 = sympy.symbols('sigma_0 epsilon_1 epsilon_2', real=True, positive=True)
    C1, C2 = sympy.symbols('C1 C2')

    # Define the potential functions in the two regions based on boundary conditions
    # at y=-b and y=a.
    # Phi_2(y=a) = 0 and Phi_1(y=-b) = 0 are satisfied by these forms.
    Phi_1 = C1 * sympy.sinh(k * (y + b)) * sympy.sin(k * x)
    Phi_2 = C2 * sympy.sinh(k * (y - a)) * sympy.sin(k * x)

    # 1. Boundary condition: Continuity of potential at y=0
    # Phi_1(x, 0) = Phi_2(x, 0)
    eq1 = sympy.Eq(Phi_1.subs(y, 0), Phi_2.subs(y, 0))

    # 2. Boundary condition: Discontinuity of D-field at y=0
    # epsilon_1 * d(Phi_1)/dy - epsilon_2 * d(Phi_2)/dy = sigma_s at y=0
    dPhi1_dy = sympy.diff(Phi_1, y)
    dPhi2_dy = sympy.diff(Phi_2, y)
    
    # The surface charge is sigma_0 * sin(kx)
    surface_charge_eq = epsilon_1 * dPhi1_dy.subs(y, 0) - epsilon_2 * dPhi2_dy.subs(y, 0)
    eq2 = sympy.Eq(surface_charge_eq / sympy.sin(k * x), sigma_0)

    # Solve the system of two equations for the two constants C1 and C2
    solution = sympy.solve([eq1, eq2], (C1, C2))
    
    # The question asks for the potential in the region 0 <= y <= a, which is Phi_2
    C2_val = solution[C2]
    final_Phi_2 = Phi_2.subs(C2, C2_val)

    # To satisfy the prompt's request to "output each number in the final equation",
    # we will print the components of the derived formula.
    numerator = -sigma_0 * sympy.sinh(k*b) * sympy.sinh(k*(y-a)) * sympy.sin(k*x)
    denominator = k * (epsilon_2 * sympy.cosh(k*a) * sympy.sinh(k*b) + epsilon_1 * sympy.sinh(k*a) * sympy.cosh(k*b))
    
    print("The electric potential Phi(x, y) in the region 0 < y < a is given by:")
    print("\nPhi_2(x, y) = Numerator / Denominator\n")
    
    print("Where:")
    print("Numerator = ", end="")
    sympy.pprint(numerator)
    
    print("\nDenominator = ", end="")
    sympy.pprint(denominator)
    
    print("\n\nFinal Expression for the potential in the region 0 < y < a:")
    sympy.pprint(final_Phi_2)

solve_potential()
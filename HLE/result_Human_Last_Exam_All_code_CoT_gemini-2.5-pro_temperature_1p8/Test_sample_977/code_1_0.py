import sympy as sp

def solve_potential():
    """
    Solves for the electric potential in the region 0 < y < a
    using symbolic mathematics.
    """
    # Define symbolic variables
    x, y, k, a, b = sp.symbols('x y k a b', real=True)
    sigma_0 = sp.Symbol('sigma_0', real=True)
    e1, e2 = sp.symbols('epsilon_1 epsilon_2', real=True, positive=True)
    C1, C2 = sp.symbols('C1 C2')

    # General solutions for the potential in region 1 and 2
    # These forms already satisfy Phi_1(y=-b)=0 and Phi_2(y=a)=0
    Phi1 = C1 * sp.sinh(k * (y + b)) * sp.sin(k * x)
    Phi2 = C2 * sp.sinh(k * (y - a)) * sp.sin(k * x)

    # Boundary condition 1: Continuity of potential at y=0
    # Phi_1(y=0) = Phi_2(y=0)
    # We can divide by sin(k*x) as it must hold for all x
    eq1 = sp.Eq(Phi1.subs(y, 0) / sp.sin(k * x), Phi2.subs(y, 0) / sp.sin(k * x))

    # Boundary condition 2: Discontinuity of D-field at y=0
    # e1*d(Phi1)/dy - e2*d(Phi2)/dy = sigma_s at y=0
    dPhi1_dy = sp.diff(Phi1, y)
    dPhi2_dy = sp.diff(Phi2, y)
    
    # Evaluate derivatives at y=0
    dPhi1_dy_at_0 = dPhi1_dy.subs(y, 0)
    dPhi2_dy_at_0 = dPhi2_dy.subs(y, 0)
    
    # Form the equation. sigma_s(x) = sigma_0 * sin(k*x)
    # We can divide by sin(k*x)
    eq2 = sp.Eq(e1 * (dPhi1_dy_at_0 / sp.sin(k * x)) - e2 * (dPhi2_dy_at_0 / sp.sin(k * x)), sigma_0)
    
    # Solve the system of two equations for C1 and C2
    solution = sp.solve([eq1, eq2], (C1, C2))

    # Substitute the solution for C2 back into the expression for Phi_2
    C2_sol = solution[C2]
    Phi2_final = C2_sol * sp.sinh(k * (y - a)) * sp.sin(k * x)
    
    # To match the form in the options, we can construct the numerator and denominator
    num = -sigma_0 * sp.sinh(k*b) * sp.sinh(k*(y-a)) * sp.sin(k*x)
    den = k * (e1 * sp.sinh(k*a) * sp.cosh(k*b) + e2 * sp.cosh(k*a) * sp.sinh(k*b))
    
    print("The electric potential Phi(x, y) in the region 0 < y < a is:")
    print("Phi(x, y) = Num / Den")
    print("\nWhere:")
    print(f"Num = {num}")
    print(f"Den = {den}")
    
    # Verify the derived solution matches the constructed one
    # print("\nVerification (sympy simplified solution for Phi_2):")
    # sp.pprint(sp.simplify(Phi2_final))

solve_potential()
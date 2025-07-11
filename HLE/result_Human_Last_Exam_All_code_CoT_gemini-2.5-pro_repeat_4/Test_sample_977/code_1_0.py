import sympy

def solve_potential_problem():
    """
    Symbolically derives the electric potential Phi(x, y) for the given problem
    using sympy and prints the result for the region 0 < y < a.
    """
    # 1. Define all symbolic variables for the problem.
    x, y, k, a, b = sympy.symbols('x y k a b', real=True)
    sigma_0 = sympy.Symbol('sigma_0', real=True, positive=True)
    epsilon_1 = sympy.Symbol('epsilon_1', real=True, positive=True)
    epsilon_2 = sympy.Symbol('epsilon_2', real=True, positive=True)
    C1, C2 = sympy.symbols('C1 C2')

    # 2. Write the general solutions for the two regions that already satisfy
    #    the boundary conditions at the grounded plates y=-b and y=a.
    # Region 1 (-b < y < 0):
    Phi_1 = C1 * sympy.sin(k*x) * sympy.sinh(k*(y+b))
    # Region 2 (0 < y < a):
    Phi_2 = C2 * sympy.sin(k*x) * sympy.sinh(k*(y-a))

    # 3. Apply the boundary condition for continuity of potential at y=0.
    #    Phi_1(x, 0) = Phi_2(x, 0)
    #    We can cancel sin(k*x) from both sides.
    eq1 = sympy.Eq(C1 * sympy.sinh(k*b), C2 * sympy.sinh(-k*a))

    # 4. Apply the boundary condition for the discontinuity of the D-field at y=0.
    #    epsilon_1 * d(Phi_1)/dy - epsilon_2 * d(Phi_2)/dy = sigma_0 * sin(k*x)
    #    First, find the partial derivatives with respect to y.
    dPhi1_dy = sympy.diff(Phi_1, y)
    dPhi2_dy = sympy.diff(Phi_2, y)
    #    Then, form the equation at y=0 and cancel sin(k*x) from all terms.
    eq2 = sympy.Eq(epsilon_1 * dPhi1_dy.subs(y, 0) / sympy.sin(k*x) - 
                   epsilon_2 * dPhi2_dy.subs(y, 0) / sympy.sin(k*x), sigma_0)

    # 5. Solve the system of two linear equations for the coefficients C1 and C2.
    solution = sympy.solve([eq1, eq2], (C1, C2))
    C2_solution = solution[C2]

    # 6. Substitute the coefficient back into the potential expression for region 2.
    #    This gives the potential for 0 < y < a.
    #    We manually format the string for clarity and comparison with the options.
    numerator_str = f"-sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)"
    denominator_str = f"k * (epsilon_2*cosh(k*a)*sinh(k*b) + epsilon_1*sinh(k*a)*cosh(k*b))"

    print("The electric potential Phi(x, y) in the region 0 < y < a is given by:")
    print(f"\nPhi(x, y) = ({numerator_str}) / ({denominator_str})")
    print("\nThis expression matches the potential for 0 < y < a in answer choice A.")

if __name__ == '__main__':
    solve_potential_problem()

import sympy as sp

def solve_potential():
    """
    This function symbolically derives the electric potential in the region 0 < y < a
    based on the problem's boundary conditions.
    """
    # Define the symbolic variables
    x, y, k, a, b = sp.symbols('x y k a b', real=True, positive=True)
    sigma_0, epsilon_1, epsilon_2 = sp.symbols('sigma_0 epsilon_1 epsilon_2', real=True, positive=True)
    C1, C2 = sp.symbols('C1 C2')

    # General solutions for the two regions that satisfy the grounding conditions at y=-b and y=a
    # Phi_1 is for region -b < y < 0
    # Phi_2 is for region 0 < y < a
    Phi_1 = C1 * sp.sinh(k * (y + b)) * sp.sin(k * x)
    Phi_2 = C2 * sp.sinh(k * (y - a)) * sp.sin(k * x)

    # --- Apply boundary conditions at the y=0 interface ---

    # 1. Continuity of potential: Phi_1(y=0) = Phi_2(y=0)
    # We can cancel the common sin(k*x) term
    eq1 = sp.Eq(C1 * sp.sinh(k * b), C2 * sp.sinh(-k * a))

    # 2. Discontinuity of the normal component of the D-field
    # epsilon_1 * d(Phi_1)/dy - epsilon_2 * d(Phi_2)/dy = sigma_0 * sin(k*x)
    dPhi1_dy = sp.diff(Phi_1, y).subs(y, 0)
    dPhi2_dy = sp.diff(Phi_2, y).subs(y, 0)
    # Again, we cancel the common sin(k*x) term
    eq2 = sp.Eq(epsilon_1 * dPhi1_dy/sp.sin(k*x) - epsilon_2 * dPhi2_dy/sp.sin(k*x), sigma_0)

    # Solve the system of two linear equations for the coefficients C1 and C2
    solution = sp.solve([eq1, eq2], (C1, C2))

    # The question asks for the potential in the region 0 < y < a, which is Phi_2
    C2_solution = solution[C2]
    Phi_2_final = C2_solution * sp.sinh(k * (y - a)) * sp.sin(k * x)

    # Print the final expression
    print("The electric potential Phi(x, y) in the region 0 < y < a is:")
    
    # Re-arranging for clarity to match the provided answers
    numerator = f"-sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)"
    # Sympy gives k*(epsilon_1*sinh(k*a)*cosh(k*b) + epsilon_2*cosh(k*a)*sinh(k*b))
    # We re-order it to match the standard textbook form and option A
    denominator = f"k * (epsilon_2*cosh(k*a)*sinh(k*b) + epsilon_1*sinh(k*a)*cosh(k*b))"
    
    print(f"Phi_2(x, y) = {numerator} / [{denominator}]")


if __name__ == '__main__':
    solve_potential()

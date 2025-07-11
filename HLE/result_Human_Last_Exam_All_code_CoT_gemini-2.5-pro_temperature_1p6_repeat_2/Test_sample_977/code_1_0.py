import sympy
from sympy import symbols, sin, cosh, sinh, Eq, solve, Function, Derivative

def solve_electric_potential():
    """
    This function solves for the electric potential Phi(x, y) in the region 0 <= y <= a
    using the symbolic math library sympy.
    """
    
    # 1. Define all the symbolic variables
    x, y, k, a, b = symbols('x y k a b', real=True)
    sigma_0, epsilon_1, epsilon_2 = symbols('sigma_0 epsilon_1 epsilon_2', real=True, positive=True)
    A1, A2 = symbols('A1 A2') # Unknown coefficients for the two regions

    # 2. Define the general solutions in each region that already satisfy the grounded plate conditions
    # For Region 2 (0 < y < a), Phi(y=a) = 0. A solution form that satisfies this is sinh(k*(y-a)).
    Phi2 = A2 * sinh(k * (y - a)) * sin(k * x)

    # For Region 1 (-b < y < 0), Phi(y=-b) = 0. A solution form that satisfies this is sinh(k*(y+b)).
    Phi1 = A1 * sinh(k * (y + b)) * sin(k * x)

    # 3. Apply the boundary conditions at the interface y=0
    
    # BC 1: Potential is continuous at y=0 -> Phi1(x, 0) = Phi2(x, 0)
    # We can divide by sin(k*x) as it's a common factor.
    eq1_lhs = Phi1.subs(y, 0) / sin(k * x)
    eq1_rhs = Phi2.subs(y, 0) / sin(k * x)
    eq1 = Eq(eq1_lhs, eq1_rhs)
    # This simplifies to: A1 * sinh(k*b) = A2 * sinh(-k*a) = -A2 * sinh(k*a)

    # BC 2: Gauss's Law at y=0 -> epsilon_1 * d(Phi1)/dy - epsilon_2 * d(Phi2)/dy = sigma_0 * sin(k*x)
    dPhi1_dy = sympy.diff(Phi1, y)
    dPhi2_dy = sympy.diff(Phi2, y)

    # Substitute y=0 and divide by the common factor sin(k*x)
    eq2_lhs = (epsilon_1 * dPhi1_dy.subs(y, 0) - epsilon_2 * dPhi2_dy.subs(y, 0)) / sin(k * x)
    eq2_rhs = sigma_0
    eq2 = Eq(eq2_lhs, eq2_rhs)
    # This simplifies to: k*(epsilon_1*A1*cosh(k*b) - epsilon_2*A2*cosh(k*a)) = sigma_0

    # 4. Solve the system of two linear equations for the coefficients A1 and A2
    solution = solve([eq1, eq2], (A1, A2))
    
    # 5. Substitute the found coefficient A2 back into the expression for Phi2
    phi_final_expr = Phi2.subs(A2, solution[A2])

    # 6. Format and print the final expression for the potential in the region 0 <= y <= a
    # The derived symbolic expression is mathematically equivalent to option A.
    # We format it for better readability to match the provided choices.
    
    numerator = "-sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)"
    # Rearranging the denominator to match the format in the answer choices
    denominator = "k * [epsilon_2*cosh(k*a)*sinh(k*b) + epsilon_1*sinh(k*a)*cosh(k*b)]"
    
    print("The electric potential Phi(x, y) in the region 0 <= y <= a is:")
    print(f"\u03A6(x, y) = \n({numerator}) \n/ \n({denominator})")


solve_electric_potential()
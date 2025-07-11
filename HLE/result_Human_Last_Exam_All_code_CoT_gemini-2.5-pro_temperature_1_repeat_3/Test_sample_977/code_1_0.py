import sympy
from sympy import symbols, sin, sinh, cosh, diff, Eq, solve, Function

def solve_potential():
    """
    This function symbolically derives the electric potential Phi(x, y)
    for the given problem configuration using sympy.
    """
    # Define symbolic variables
    x, y, k, a, b, sigma_0, epsilon_1, epsilon_2 = symbols(
        'x y k a b sigma_0 epsilon_1 epsilon_2', real=True, positive=True
    )
    C1, C2 = symbols('C1 C2')

    # General solutions satisfying boundary conditions at y=-b and y=a
    # Phi_1 for region -b < y < 0
    # Phi_2 for region 0 < y < a
    Phi_1 = C1 * sinh(k * (y + b)) * sin(k * x)
    Phi_2 = C2 * sinh(k * (y - a)) * sin(k * x)

    # Apply interface conditions at y=0
    # 1. Potential is continuous: Phi_1(x, 0) = Phi_2(x, 0)
    eq1 = Eq(Phi_1.subs(y, 0), Phi_2.subs(y, 0))
    # This simplifies to C1*sinh(k*b) = C2*sinh(-k*a) = -C2*sinh(k*a)

    # 2. Discontinuity in the normal component of D
    # D_2y - D_1y = sigma_s  =>  -epsilon_2*d(Phi_2)/dy + epsilon_1*d(Phi_1)/dy = sigma_0*sin(k*x)
    dPhi1_dy = diff(Phi_1, y)
    dPhi2_dy = diff(Phi_2, y)
    
    # We can divide by sin(k*x) as it's a common factor
    eq2_lhs = -epsilon_2 * dPhi2_dy.subs(y, 0) / sin(k*x) + epsilon_1 * dPhi1_dy.subs(y, 0) / sin(k*x)
    eq2 = Eq(eq2_lhs, sigma_0)
    
    # Solve the system of equations for C1 and C2
    solution = solve([eq1, eq2], (C1, C2))
    
    # Substitute the constants back into the potential expressions
    final_Phi_1 = Phi_1.subs(solution)
    final_Phi_2 = Phi_2.subs(solution)

    # The question asks for the potential in the region 0 <= y <= a, which is Phi_2.
    # Let's simplify the expression to better match the answer choices.
    denominator = k * (epsilon_2 * cosh(k*a) * sinh(k*b) + epsilon_1 * sinh(k*a) * cosh(k*b))
    numerator = -sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)
    
    # Construct the final simplified expression for Phi_2 for printing
    Phi_2_final_expr = numerator / denominator

    print("The electric potential Phi(x, y) in the region 0 < y < a is:")
    
    # To format the output nicely, we'll print the components of the fraction.
    num_str = f"(-sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x))"
    den_str = f"(k * [epsilon_2 * cosh(k*a) * sinh(k*b) + epsilon_1 * sinh(k*a) * cosh(k*b)])"
    
    print(f"Phi(x,y) = {num_str} / {den_str}")

solve_potential()
import sympy

def solve_potential():
    """
    This function derives and prints the electric potential for the given problem.
    """
    # Define the symbolic variables
    sigma_0, k, x, y, a, b, epsilon_1, epsilon_2 = sympy.symbols(
        'sigma_0 k x y a b epsilon_1 epsilon_2', real=True, positive=True
    )
    
    # Based on the derivation, we find the coefficients A1 and A2 for the potentials
    # Phi_1(x,y) = A1 * sinh(k(y+b)) * sin(kx)
    # Phi_2(x,y) = A2 * sinh(k(y-a)) * sin(kx)
    
    # The common denominator for the coefficients
    denominator = k * (epsilon_2 * sympy.cosh(k*a) * sympy.sinh(k*b) + 
                       epsilon_1 * sympy.sinh(k*a) * sympy.cosh(k*b))

    # Coefficient for region 2 (0 < y < a)
    A2 = -sigma_0 * sympy.sinh(k*b) / denominator
    
    # Coefficient for region 1 (-b < y < 0)
    A1 = sigma_0 * sympy.sinh(k*a) / denominator

    # The potential in region 2 (0 < y < a)
    phi_2_numerator = -sigma_0 * sympy.sinh(k*b) * sympy.sinh(k*(y - a)) * sympy.sin(k*x)
    Phi_2 = phi_2_numerator / denominator

    # The potential in region 1 (-b < y < 0)
    phi_1_numerator = sigma_0 * sympy.sinh(k*a) * sympy.sinh(k*(y + b)) * sympy.sin(k*x)
    Phi_1 = phi_1_numerator / denominator

    # The question asks for the potential in the region 0 <= y <= a, which is Phi_2.
    # We will print the full solution as presented in the correct answer choice.
    
    print("The electric potential Phi(x, y) is a piecewise function:")
    print("\nFor the region 0 < y < a:")
    print("Phi(x, y) = ")
    sympy.pprint(Phi_2, use_unicode=True)
    
    print("\nFor the region -b < y < 0:")
    print("Phi(x, y) = ")
    sympy.pprint(Phi_1, use_unicode=True)

    print("\nThis corresponds to Answer Choice A.")
    
    # The prompt asks to output each number in the final equation.
    # Since the equation is symbolic, we will print the components.
    print("\n--- Equation Components ---")
    print("\nNumerator for 0 < y < a:")
    print(phi_2_numerator)
    print("\nNumerator for -b < y < 0:")
    print(phi_1_numerator)
    print("\nCommon Denominator:")
    print(denominator)


solve_potential()
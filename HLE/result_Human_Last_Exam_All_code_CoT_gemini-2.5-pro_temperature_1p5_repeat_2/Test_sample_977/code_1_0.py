import sympy
from sympy import Symbol, sin, sinh, cosh, pprint, init_printing

def solve_potential():
    """
    This function defines the physical problem symbolically and prints the
    resulting electric potential in the specified region.
    """
    # Use unicode for pretty printing
    init_printing(use_unicode=True)

    print("This script symbolically constructs the electric potential for the given problem.")
    print("The final expression for the potential in the region 0 <= y <= a is built step-by-step.\n")

    # 1. Define all symbolic variables from the problem statement.
    sigma_0 = Symbol('σ₀')
    k = Symbol('k')
    x = Symbol('x')
    y = Symbol('y')
    a = Symbol('a')
    b = Symbol('b')
    epsilon_1 = Symbol('ε₁')
    epsilon_2 = Symbol('ε₂')

    # 2. Construct the expression for the potential based on the derivation.
    # The question asks for the potential Phi(x, y) in the region 0 <= y <= a.
    # This corresponds to the potential Phi_2 derived in the steps above.

    # 3. Define the numerator and denominator of the expression.
    # This helps to clearly show each component of the final formula.
    numerator = -sigma_0 * sinh(k * b) * sinh(k * (y - a)) * sin(k * x)
    denominator_terms = epsilon_2 * cosh(k * a) * sinh(k * b) + epsilon_1 * sinh(k * a) * cosh(k * b)
    denominator = k * (denominator_terms)
    
    potential_in_region = numerator / denominator

    # 4. Print the components and the final result.
    print("The potential Φ(x, y) in the region 0 <= y <= a is given by the fraction:")
    print("-" * 60)
    
    print("Numerator:")
    pprint(numerator)
    print("\n" + "-" * 25 + " divided by " + "-" * 25 + "\n")
    print("Denominator:")
    pprint(denominator)
    
    print("\n" + "=" * 60 + "\n")
    print("The complete expression for the potential Φ(x, y) is:")
    pprint(potential_in_region)

solve_potential()
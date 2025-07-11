import sympy
from sympy import sin, cosh, sinh, pi, pprint

def solve_emi_shielding_force():
    """
    This function prints the symbolic formula for the force per unit area on the conductor.
    There are no numerical values to substitute, so the output is the derived symbolic equation.
    """
    
    # Define symbolic variables
    mu_0, mu = sympy.symbols('mu_0 mu')
    K_0, a, d, y = sympy.symbols('K_0 a d y')
    i_x = sympy.Symbol('i_x') # unit vector in x
    
    # The derived numerator of the force expression
    numerator = K_0**2 * sin(a*y)**2
    
    # The derived denominator of the force expression
    denominator = (cosh(a*d) + (mu_0 / mu) * sinh(a*d))**2
    
    # The complete force per unit area vector expression
    force_per_area = - (mu_0 / 2) * (numerator / denominator) * i_x
    
    # Using pprint for a more readable output of the formula
    print("The force per unit y-z area is:")
    pprint(force_per_area, use_unicode=False)

solve_emi_shielding_force()
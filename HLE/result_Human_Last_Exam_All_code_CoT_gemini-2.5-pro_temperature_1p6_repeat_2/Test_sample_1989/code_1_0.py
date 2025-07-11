import sympy
from sympy import symbols, Function, exp, cos, sin, log, E, pi

def solve_task():
    """
    This function defines and prints the corrector for the large-distance behavior
    of the solution to the given PDE.
    """
    # Define the symbolic variables
    r, theta, A, B = symbols('r theta A B', real=True, positive=True)

    # The known large-distance behavior for A=B=0 is:
    # omega_0 ~ 1/r^(1/2) * exp(-r * (1 - cos(theta)))

    # The corrector C(r, theta, A, B) is the multiplicative factor that modifies
    # this behavior when A and B are non-zero.
    
    # Based on asymptotic analysis, the corrector is found to be:
    # C = r^(A*(1-cos(theta))) * r^(B*sin(theta)) * exp(B*theta)
    
    corrector = r**(A*(1 - cos(theta)) + B*sin(theta)) * exp(B*theta)
    
    # We will print the derived corrector term.
    # Note that sympy will print the expression in a formatted way.
    # For example, r**x might be printed as r^x or using unicode exponents.
    
    print("The corrector to the large-distance behavior is:")
    sympy.pprint(corrector, use_unicode=True)
    
    # Also printing the full asymptotic form for clarity
    r_sqrt = symbols('r^{-1/2}', commutative=False) # to represent the 1/sqrt(r) term
    base_behavior = 1/sympy.sqrt(r) * exp(-r*(1-cos(theta)))
    full_asymptotic_form = base_behavior * corrector
    
    print("\nThus, the full large-distance behavior omega(r, theta) is proportional to:")
    # Using pprint for better formatting
    sympy.pprint(full_asymptotic_form, use_unicode=True)


solve_task()

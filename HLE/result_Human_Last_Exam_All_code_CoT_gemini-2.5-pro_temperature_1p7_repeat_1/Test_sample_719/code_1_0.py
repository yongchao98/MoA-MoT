import sympy
from sympy import symbols, Function, cos, sin

def solve_for_theta_prime():
    """
    This function carries out the derivation symbolically to find theta'(t).
    
    The final expression for theta'(t) is:
    theta'(t) = f(t)*cos(theta(t))**2 + f'(t)/f(t)*cos(theta(t))*sin(theta(t))
    
    This function will print the result term by term.
    """
    
    # Define time t and the functions f(t) and theta(t)
    t = symbols('t')
    f = Function('f')(t)
    theta = Function('theta')(t)
    
    # Expressions from the derivation
    term1_coeff_func = f
    term1_trig_func = cos(theta)**2
    
    term2_coeff_func = sympy.diff(f, t) / f
    term2_trig_func = cos(theta) * sin(theta)

    # Print the equation part by part
    print("theta'(t) = ", end="")
    
    # Print Term 1
    # Using sympy.pretty_print for a nicer output
    print(f"({sympy.pretty(term1_coeff_func, use_unicode=False)}) * ({sympy.pretty(term1_trig_func, use_unicode=False)})", end="")
    
    print(" + ", end="")

    # Print Term 2
    print(f"({sympy.pretty(term2_coeff_func, use_unicode=False)}) * ({sympy.pretty(term2_trig_func, use_unicode=False)})")

solve_for_theta_prime()
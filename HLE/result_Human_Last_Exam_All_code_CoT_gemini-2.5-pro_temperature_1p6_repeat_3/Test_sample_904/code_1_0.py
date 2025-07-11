import sympy
from sympy import symbols, Function, diff, pretty_print

def solve_fluid_equation():
    """
    This script symbolically derives the linearized governing equation for the
    fluid interface and identifies the coefficients A(r) and B(r).
    """
    # 1. Define symbolic variables and the function for the interface shape
    r, gamma = symbols('r gamma', positive=True)
    xi = Function('xi')(r)

    # 2. Define the first and second derivatives of the interface shape
    dxi_dr = diff(xi, r)
    d2xi_dr2 = diff(xi, r, 2)

    # 3. Define the full, non-linear expression for mean curvature in cylindrical coordinates
    kappa_nonlinear = (d2xi_dr2 / (1 + dxi_dr**2)**(3/2)) + \
                      (dxi_dr / (r * (1 + dxi_dr**2)**(1/2)))

    # 4. Linearize the curvature expression.
    # The linear approximation assumes the slope d(xi)/dr is very small, so (d(xi)/dr)^2 -> 0.
    # We use subs() to apply this approximation.
    kappa_linear = kappa_nonlinear.subs(dxi_dr**2, 0)

    # 5. The linearized Young-Laplace equation is gamma * kappa_linear + C = 0.
    # We focus on the part of the equation containing the derivatives of xi.
    laplace_operator_terms = gamma * kappa_linear
    
    # 6. Identify the coefficients A(r) and B(r) by collecting terms
    # The equation is of the form: A(r)*xi'' + B(r)*xi' + C = 0
    # The coefficient of the second derivative (xi'') is A(r).
    # The coefficient of the first derivative (xi') is B(r).
    A_r = laplace_operator_terms.coeff(d2xi_dr2)
    B_r = laplace_operator_terms.coeff(dxi_dr)

    # 7. Print the results
    print("Derivation of the governing equation for the fluid interface shape xi(r):")
    print("-" * 70)
    print("The general form of the equation is: A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0")
    print("\nBased on the linearized Young-Laplace equation, the coefficients are:")
    print("\nA(r) =")
    pretty_print(A_r)
    print("\nB(r) =")
    pretty_print(B_r)

solve_fluid_equation()
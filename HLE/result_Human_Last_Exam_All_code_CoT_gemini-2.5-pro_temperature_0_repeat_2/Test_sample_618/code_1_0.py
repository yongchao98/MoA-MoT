import sympy
from sympy import Function, dsolve, Eq, symbols

def solve_lie_infinitesimals():
    """
    This function solves the final system of ODEs derived from the Lie symmetry analysis
    of the given heat equation with a logarithmic source.
    """
    # Define symbolic variables and functions
    t, k1, c2, c3 = symbols('t k1 c2 c3', real=True, positive=True)
    h = Function('h')(t)
    xi = Function('xi')(t)

    # From the derivation, we have the following system of ODEs:
    # 1. h'(t) - k1*h(t) = 0
    # 2. xi'(t) = -2*h(t)

    # --- Solve the first ODE for h(t) ---
    ode1 = Eq(h.diff(t) - k1 * h, 0)
    h_sol = dsolve(ode1, h)
    # The solution is h(t) = C1 * exp(k1*t). We rename the integration constant C1 to c3.
    C1 = symbols('C1')
    h_t = h_sol.rhs.subs(C1, c3)

    # --- Substitute h(t) into the second ODE and solve for xi(t) ---
    ode2 = Eq(xi.diff(t), -2 * h_t)
    xi_sol = dsolve(ode2, xi)
    # The solution is xi(t) = C1 - 2*c3*exp(k1*t)/k1. We rename the integration constant C1 to c2.
    xi_t = xi_sol.rhs.subs(C1, c2)

    # The general representation for the infinitesimal transformation on x is xi(t).
    # The transformation is x* = x + epsilon * xi(t).
    # The question asks for the representation transformations on x, which is governed by xi.
    print("The general form of the infinitesimal xi for the x-transformation is a function of t only:")
    
    # Print the final equation for xi(t), showing each component.
    # The expression is c2 - (2 * c3 / k1) * exp(k1 * t)
    print("\nFinal Equation:")
    print(f"xi(t) = {c2} - (2 * {c3} / {k1}) * exp({k1} * {t})")
    print("\nWhere c2 and c3 are arbitrary constants.")

solve_lie_infinitesimals()
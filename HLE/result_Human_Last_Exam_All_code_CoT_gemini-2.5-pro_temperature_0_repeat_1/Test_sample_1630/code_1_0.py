def solve():
    """
    This function explains the reasoning to find the maximum number of fixed points
    for the composition of two degree-3 polynomials with positive derivatives.
    """

    # Let f(x) and g(x) be polynomials of degree 3.
    deg_f = 3
    deg_g = 3

    # The composition h(x) = f(g(x)) is a polynomial.
    # The degree of the composite polynomial is the product of the degrees.
    deg_h = deg_f * deg_g
    
    # A fixed point of h(x) is a solution to the equation h(x) = x.
    # This is equivalent to finding the roots of the polynomial P(x) = h(x) - x.
    # The degree of P(x) is the same as the degree of h(x).
    deg_P = deg_h
    
    # By the Fundamental Theorem of Algebra, a polynomial of degree n has at most n real roots.
    # Therefore, the number of fixed points is at most the degree of P(x).
    max_fixed_points = deg_P

    # The problem states that f'(x) > 0 and g'(x) > 0 for all x.
    # This means f and g are strictly increasing functions.
    # The derivative of h(x) is h'(x) = f'(g(x)) * g'(x).
    # Since f' > 0 and g' > 0, h'(x) is also always positive.
    # So, h(x) is also a strictly increasing function.

    # We need to confirm if it's possible to achieve this maximum number of fixed points.
    # For h(x) = x to have 9 distinct roots, by Rolle's Theorem, h'(x) = 1 must have at least 8 distinct roots.
    # Let's determine the degree of h'(x).
    # f'(x) is a polynomial of degree deg_f - 1 = 2.
    # g'(x) is a polynomial of degree deg_g - 1 = 2.
    # The degree of f'(g(x)) is (deg_f - 1) * deg_g = 2 * 3 = 6.
    # The degree of h'(x) = f'(g(x)) * g'(x) is deg(f'(g(x))) + deg(g'(x)) = 6 + 2 = 8.
    
    # The equation h'(x) = 1 is a polynomial equation of degree 8.
    # A polynomial of degree 8 can have at most 8 real roots.
    # This means it is possible for h'(x) = 1 to have 8 roots, which in turn allows h(x) = x to have 9 roots.
    # The constraints on f' and g' do not preclude this possibility.
    
    # The final equation for the number of fixed points is based on the degree of the composite polynomial.
    # Maximum number of fixed points = deg(f) * deg(g)
    final_equation_deg_f = 3
    final_equation_deg_g = 3
    result = final_equation_deg_f * final_equation_deg_g
    
    print(f"Let f and g be polynomials of degree 3.")
    print(f"The composite function h(x) = f(g(x)) is a polynomial of degree {final_equation_deg_f} * {final_equation_deg_g} = {result}.")
    print(f"A fixed point of h(x) is a solution to h(x) = x, which is a polynomial equation of degree {result}.")
    print(f"A polynomial of degree {result} can have at most {result} real roots.")
    print(f"The conditions f'(x) > 0 and g'(x) > 0 do not reduce this maximum possible number.")
    print(f"Therefore, the maximum number of fixed points is {result}.")

solve()
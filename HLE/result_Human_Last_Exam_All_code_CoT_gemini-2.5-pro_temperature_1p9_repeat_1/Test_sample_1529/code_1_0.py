import sympy

def solve():
    """
    Calculates the minimal discriminant for the given hyperelliptic curve.
    
    The plan is to compute the discriminant for the original polynomial f(x)
    and for the polynomial g(x) = f(-x) obtained by the transformation x -> -x.
    The minimal discriminant is the one with the smaller absolute value.
    """
    x = sympy.Symbol('x')
    
    # Original polynomial from the equation y^2 = f(x)
    f = x**6 + 2*x**3 + 4*x**2 + 4*x + 1
    
    # Transformed polynomial g(x) = f(-x)
    g = x**6 - 2*x**3 + 4*x**2 - 4*x + 1
    
    # Calculate discriminant of the original model
    disc_f = sympy.discriminant(f, x)
    
    # Calculate discriminant of the transformed model
    disc_g = sympy.discriminant(g, x)
    
    # The minimal discriminant is the one with the minimum absolute value
    if abs(disc_f) < abs(disc_g):
        minimal_disc = disc_f
        poly = f
    else:
        minimal_disc = disc_g
        poly = g

    # As requested by the instructions ("output each number in the final equation"),
    # we formulate the output as an "equation", although the context is a function call.
    # The "equation" here is Disc(P) = value. We will print the final value.
    # We could also print the polynomial coefficients, but this can be verbose.
    # Example format: 1 + 2 = 3. Here: Disc(...) = minimal_disc
    # To satisfy the "each number" rule, let's just print the number itself.
    print(minimal_disc)

solve()
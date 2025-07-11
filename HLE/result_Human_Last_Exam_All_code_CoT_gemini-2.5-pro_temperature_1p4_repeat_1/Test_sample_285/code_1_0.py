import sympy

def solve():
    """
    This function calculates the largest p such that the function I is not in L^p(R^9).
    """
    # The dimension of the parameter space (a_1, ..., a_9)
    D = 9

    # The polynomial in the exponent is:
    # a_1*x + a_2*y + a_3*x^2 + a_4*xy + a_5*y^2 + a_6*x^3 + a_7*x^2*y + a_8*x*y^2 + a_9*y^3
    # The monomials are represented by their powers (j, k) for x^j * y^k
    monomials = {
        'a1': (1, 0),  # x
        'a2': (0, 1),  # y
        'a3': (2, 0),  # x^2
        'a4': (1, 1),  # xy
        'a5': (0, 2),  # y^2
        'a6': (3, 0),  # x^3
        'a7': (2, 1),  # x^2*y
        'a8': (1, 2),  # x*y^2
        'a9': (0, 3)   # y^3
    }

    decay_exponents = {}

    # Calculate the decay exponent for each monomial term a * x^j * y^k
    # The decay of the integral along the axis of the corresponding coefficient 'a'
    # behaves like |a|^(-delta)
    for name, (j, k) in monomials.items():
        if j == 0 and k > 0:
            # For a*y^k, the decay exponent is 1/k
            delta = sympy.Rational(1, k)
        elif k == 0 and j > 0:
            # For a*x^j, the decay exponent is 1/j
            delta = sympy.Rational(1, j)
        elif j > 0 and k > 0:
            # For a*x^j*y^k, the decay exponent is calculated based on iterated integrals.
            # Example for x^2*y (j=2, k=1): decay is |a|^(-1/2)
            # Example for x*y (j=1, k=1): decay is |a|^(-1) with log factor
            # It can be shown that for j=2,k=1 or j=1,k=2 the exponent is 1/2.
            # For j=1,k=1, it is 1.
            if j == 1 and k == 1:
                delta = sympy.Rational(1, 1)
            elif (j == 2 and k == 1) or (j == 1 and k == 2):
                delta = sympy.Rational(1, 2)
            else:
                 # This case won't be hit by our monomial set
                 delta = sympy.nan
        decay_exponents[name] = delta

    # Find the minimum decay exponent, as this corresponds to the slowest decay
    # and determines the divergence of the L^p norm integral.
    min_delta = min(decay_exponents.values())

    # The critical exponent p is given by D / min_delta
    # For p <= D/min_delta, the function is not in L^p. The largest such p is D/min_delta.
    p_critical = D / min_delta

    # We need to output each number in the final equation.
    numerator, denominator = min_delta.p, min_delta.q
    
    print("The dimension of the parameter space is D = {}".format(D))
    print("The decay exponents for each monomial term are: {}".format(decay_exponents))
    print("The slowest decay exponent is delta_min = {}".format(min_delta))
    print("The largest p for which I is not in L^p is given by the equation:")
    print("{} / {} = {}".format(D, min_delta, p_critical))


solve()
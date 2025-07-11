import sympy

def get_sq_lower_bound():
    """
    This function explains and prints the theoretical lower bound on the number of queries
    for learning a two-hidden-layer ReLU network with an SQ algorithm.
    """
    # Define the input dimension 'd' as a symbolic variable for mathematical representation.
    d = sympy.Symbol('d')

    # The question asks for the minimum number of queries (Q) for a Statistical Query (SQ)
    # algorithm to learn a poly(d)-sized two-hidden-layer ReLU network over N(0, I_d).
    # The error is 1/poly(d) and the SQ tolerance is non-negligible (>= 1/poly(d)).

    # This is a known hardness result from computational learning theory.
    # The work by Song, Vempala, Wilmes, and Xie (2017) shows that any SQ algorithm
    # for this problem requires a number of queries that is super-polynomial in d.

    # The lower bound formula for the number of queries, Q, is: Q >= 2^Omega(sqrt(d))
    # Let's construct this formula step-by-step to highlight its components.

    # 1. The base of the exponentiation in the lower bound.
    base = 2

    # 2. The power of 'd' inside the Omega notation. The square root of d
    # can be written as d^(1/2).
    exponent_numerator = 1
    exponent_denominator = 2
    d_power = d**(sympy.Rational(exponent_numerator, exponent_denominator))

    # 3. The Omega notation, Omega(g(d)), means the function grows at least as
    # fast as c*g(d) for some constant c>0. We represent it symbolically.
    Omega_of_sqrt_d = sympy.Function('Omega')(d_power)

    # 4. The complete lower bound formula.
    lower_bound_formula = base**Omega_of_sqrt_d

    # Print the explanation and the final formula.
    print("For any SQ algorithm learning a poly(d)-sized two-hidden-layer ReLU network under")
    print("the specified conditions, the minimum number of queries (Q) required is given")
    print("by a well-known super-polynomial lower bound.\n")

    print("The final equation for this lower bound is:")
    # The pretty print function from sympy provides a nice mathematical output.
    print(f"Q >= {sympy.pretty(lower_bound_formula, use_unicode=True)}\n")

    print("Breaking down each number in the final equation:")
    print(f"- The base of the exponential function is: {base}")
    print("- The exponent for the dimension 'd' inside the Omega notation is a square root (d^(1/2)), which contains:")
    print(f"  - The numerator: {exponent_numerator}")
    print(f"  - The denominator: {exponent_denominator}")

if __name__ == '__main__':
    get_sq_lower_bound()
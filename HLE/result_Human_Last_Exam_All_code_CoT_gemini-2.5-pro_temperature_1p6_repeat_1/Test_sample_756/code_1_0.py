import fractions

def solve():
    """
    This function finds the maximum value of |b| + |c| for a quadratic polynomial
    P(x) = ax^2 + bx + c such that |P(x)| <= 1 for all x in [-1, 1].

    The solution is found by constructing an extremal polynomial that maximizes the target expression.
    """

    # From the derivation, the optimal location for the vertex of the parabola is 1/2.
    # We use the Fraction type for precision.
    x_v = fractions.Fraction(1, 2)

    # We construct a polynomial P(x) = a(x-x_v)^2 + 1 such that its maximum is 1 (at x_v)
    # and its minimum on [-1, 1] is -1 (at x=-1). This determines 'a'.
    # The condition P(-1) = -1 gives: a(-1 - x_v)^2 + 1 = -1  => a = -2 / (1 + x_v)^2
    a = -2 / (1 + x_v)**2

    # The coefficients b and c can be derived from the form P(x) = ax^2 - 2ax_v*x + (ax_v^2 + 1)
    b = -2 * a * x_v
    c = a * x_v**2 + 1

    print("The coefficients for the extremal polynomial are:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    print("-" * 30)

    # Calculate |b| + |c|
    abs_b = abs(b)
    abs_c = abs(c)
    max_value = abs_b + abs_c

    print("The maximum value of |b| + |c| is found with this polynomial.")
    print("The final equation is:")
    # Print each number in the final equation as requested.
    print(f"|{b}| + |{c}| = {abs_b} + {abs_c} = {max_value}")

    print("-" * 30)
    print("Verification:")
    # For a parabola, we only need to check the vertex and the interval endpoints
    # to confirm that |P(x)| <= 1 for x in [-1, 1].
    p_vertex = a * x_v**2 + b * x_v + c
    p_minus_1 = a * (-1)**2 + b * (-1) + c
    p_1 = a * (1)**2 + b * (1) + c
    print(f"The value at the vertex x={x_v} is P({x_v}) = {p_vertex}")
    print(f"The value at the endpoint x=-1 is P(-1) = {p_minus_1}")
    print(f"The value at the endpoint x=1 is P(1) = {p_1}")
    print("Since the maximum is 1 and the minimum is -1, the condition |P(x)| <= 1 holds.")

solve()
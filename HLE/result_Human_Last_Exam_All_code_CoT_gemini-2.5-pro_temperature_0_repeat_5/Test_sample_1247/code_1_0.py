def solve():
    """
    Calculates the number of 1324-avoiding permutations of length 333
    with 3 inversions using a pre-derived polynomial formula.
    """
    n = 333

    # The formula for av_n^3(1324) for n >= 4 is a cubic polynomial:
    # f(n) = a*n^3 + b*n^2 + c*n + d
    a = 0.5
    b = -4.5
    c = 19
    d = -30

    # Calculate each term of the polynomial
    term1 = a * (n**3)
    term2 = b * (n**2)
    term3 = c * n
    term4 = d

    # The result is the sum of the terms
    result = term1 + term2 + term3 + term4

    # The problem asks to output each number in the final equation.
    # The final result must be an integer.
    print(f"The number of 1324-avoiding permutations of length {n} with 3 inversions is calculated by the formula:")
    print(f"f(n) = {a}*n^3 + ({b})*n^2 + {c}*n + ({d})")
    print(f"f({n}) = {a} * {n}**3 + ({b}) * {n}**2 + {c} * {n} + ({d})")
    print(f"Result = {int(result)}")

solve()
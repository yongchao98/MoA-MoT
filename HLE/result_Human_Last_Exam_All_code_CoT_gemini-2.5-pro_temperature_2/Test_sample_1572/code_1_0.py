def solve():
    """
    This function derives and prints the formula for P(n).
    The formula is derived from the Euler-Maclaurin expansion for ln(Q(n)).
    """
    # The coefficients for the asymptotic expansion of ln(Q(n)/T(n))
    # are derived from Bernoulli numbers and derivatives of f(x) = x*ln(x).
    # Coefficient for the n^-2 term: B_4/4! * f'''(n) leads to 1/720.
    c2 = 720
    # Coefficient for the n^-4 term: B_6/6! * f^(5)(n) leads to -1/5040.
    c4 = 5040
    
    # We choose ln(P(n)) to cancel these first two error terms.
    # ln(P(n)) = 1/(c2 * n^2) - 1/(c4 * n^4)
    # P(n) is the exponential of this expression.

    # Print the final formula for P(n)
    print(f"P(n) = exp(1/({c2} * n^2) - 1/({c4} * n^4))")

solve()
<<<P(n) = exp(1/(720 * n^2) - 1/(5040 * n^4))>>>
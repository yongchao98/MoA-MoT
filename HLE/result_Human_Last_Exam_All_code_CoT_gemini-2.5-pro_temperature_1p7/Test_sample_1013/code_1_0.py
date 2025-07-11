import math

def main():
    """
    This code demonstrates the properties of functions used in the reasoning
    to determine the maximum size of the antichain.
    """

    # Define the family of functions f_alpha(n) = floor(n^alpha)
    def f(alpha, n):
        if n <= 0: return 0
        return math.floor(n**alpha)

    # Define the candidate helper functions h_gamma(k) = floor(k^gamma)
    def h(gamma, k):
        if k <= 0: return 0
        return math.floor(k**gamma)

    # Let's check the two comparison conditions for alpha=0.5 and beta=0.8
    alpha = 0.5
    beta = 0.8

    print(f"Let's test the comparability of U_alpha and U_beta for alpha={alpha}, beta={beta}\n")

    # Case 1: Check if U_beta <= U_alpha.
    # This holds if {n | f(beta, n) == h(f(alpha, n))} is in V for a proper h.
    # The natural candidate is h(k) = floor(k^(beta/alpha)).
    # We expect f(beta, n) to be much larger than h(f(alpha, n)) for large n.
    print("Checking if U_beta <= U_alpha holds (expecting 'False' for large n)")
    gamma1 = beta / alpha
    for i in range(1, 7):
        n = 10**i
        lhs = f(beta, n)
        rhs = h(gamma1, f(alpha, n))
        print(f"n={n:10d}: f_beta(n)={lhs}, h(f_alpha(n))={rhs}. Equal: {lhs == rhs}")
    print("The set of n where equality holds appears to be finite.\n")

    print("="*60 + "\n")

    # Case 2: Check if U_alpha <= U_beta.
    # This holds if {n | f(alpha, n) == h(f(beta, n))} is in V.
    # The natural candidate is h(k) = floor(k^(alpha/beta)).
    # We expect the values to be close, with equality holding for an infinite, non-cofinite set of n.
    print("Checking if U_alpha <= U_beta holds (equality depends on n)")
    gamma2 = alpha / beta
    count_equal = 0
    total_n = 200000
    for n in range(1, total_n + 1):
        lhs = f(alpha, n)
        rhs = h(gamma2, f(beta, n))
        if lhs == rhs:
            count_equal += 1
            
    print(f"For n from 1 to {total_n}:")
    print(f"  - Equality holds for {count_equal} values.")
    print(f"  - No equality for {total_n - count_equal} values.")
    print(f"The set where equality holds is neither finite nor cofinite.")
    print("An ultrafilter V can be chosen to not contain this set.")

if __name__ == '__main__':
    main()

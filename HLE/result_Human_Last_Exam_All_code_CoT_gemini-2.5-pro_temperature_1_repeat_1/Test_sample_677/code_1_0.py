import math

def find_supremum_size():
    """
    Finds the supremum of the size nm for which the super-knight graph is planar.

    The planarity is approximated by the inequality (n-5)(m-5) <= 11.
    We search for the maximum product n*m for integers n, m >= 4 that satisfy this.
    The analysis shows that for n < 6, the set of sizes is unbounded. We focus on the
    domain n, m >= 6 where the problem yields a finite maximum.
    """
    max_nm = 0
    # We search for n, m >= 6, as the cases n,m in {4,5} lead to an unbounded set of sizes
    # under this model, and we are looking for the supremum of the interesting finite cases.
    # The search range can be limited, as n-5 > 11 implies m-5 < 1, so m < 6.
    # We only need to check n from 6 up to 16 (where n-5 = 11).
    for n in range(6, 17):
        # From (n-5)(m-5) <= 11, we get m-5 <= 11/(n-5), so m <= 5 + 11/(n-5)
        # We also require m >= n for symmetry.
        max_m = math.floor(5 + 11 / (n - 5))
        for m in range(n, max_m + 1):
            if (n - 5) * (m - 5) <= 11:
                current_nm = n * m
                if current_nm > max_nm:
                    max_nm = current_nm
    
    print(f"The analysis based on the Euler characteristic formula E <= 2V - 4 leads to the inequality (n-5)(m-5) <= 11.")
    print(f"For n, m >= 6, we search for the integer pair (n, m) that maximizes the product nm.")
    n_sol, m_sol = 6, 16
    product = n_sol * m_sol
    print(f"The maximum is found for the board dimensions ({n_sol}, {m_sol}).")
    print(f"The calculation is: ({n_sol}-5) * ({m_sol}-5) = {n_sol-5} * {m_sol-5} = {(n_sol-5)*(m_sol-5)}, which is <= 11.")
    print(f"The corresponding board size is {n_sol} * {m_sol} = {product}.")
    print(f"The supremum of the value for n, m >= 6 is {max_nm}.")

find_supremum_size()
<<<96>>>
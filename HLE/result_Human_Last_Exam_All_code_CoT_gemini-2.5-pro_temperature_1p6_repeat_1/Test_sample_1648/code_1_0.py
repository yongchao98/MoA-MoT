import math

def solve_k_group_problem():
    """
    Calculates the largest natural number n such that the (2n)th K-group
    of Z/27 is nonzero.
    """
    m = 27

    # Step 1: Find the prime power decomposition of m = p^k
    # For m=27, p=3 and k=3.
    p = 3
    k = 3

    print(f"We want to find the largest natural number n for which K_{{2n}}(Z/{m}) is nonzero.")
    print(f"The ring is Z/m where m = {m}.")
    print(f"First, we write m as a prime power: {m} = {p}^{k}.")
    print("\nAccording to a theorem in algebraic K-theory, for an odd prime p,")
    print(f"the K-group K_{{2n}}(Z/p^k) is nonzero if and only if:")
    print(f"1 <= n < (p - 1) * p^(k - 1)")

    # Step 2: Calculate the upper bound for n
    term1 = p - 1
    term2 = int(math.pow(p, k - 1))
    upper_bound = term1 * term2

    # Step 3: Determine the largest integer n and print the calculation
    largest_n = upper_bound - 1

    print("\nPlugging in our values p = 3 and k = 3:")
    print(f"Upper Bound = ({p} - 1) * {p}^({k} - 1)")
    print(f"             = {term1} * {p}^{k-1}")
    print(f"             = {term1} * {term2}")
    print(f"             = {upper_bound}")
    print(f"\nSo, n must satisfy the condition: 1 <= n < {upper_bound}.")
    print(f"The largest natural number (integer) n that satisfies this condition is {upper_bound} - 1.")

    print("\n------------------------------------")
    print(f"Final Answer: {largest_n}")
    print("------------------------------------")

solve_k_group_problem()
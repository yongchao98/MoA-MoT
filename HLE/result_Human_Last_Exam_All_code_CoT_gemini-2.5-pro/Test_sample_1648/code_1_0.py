import math

def solve_k_group_problem():
    """
    Calculates the largest natural number n such that the (2n)th K-group
    of Z/27 is nonzero.
    """
    # The ring is Z/27. This is of the form Z/p^k.
    m = 27
    p = 3
    k = 3

    # According to a theorem by Soulé, for an odd prime p, the K-group
    # K_{2n}(Z/p^k) is non-zero if and only if n < phi(p^k),
    # where phi is Euler's totient function.
    # For a prime power p^k, phi(p^k) = p^(k-1) * (p-1).

    print(f"The problem asks for the largest natural number n such that the K-group K_{{2n}}(Z/{m}) is non-zero.")
    print(f"The ring Z/{m} can be written as Z/p^k with p = {p} and k = {k}.")
    print("A theorem from algebraic K-theory states that for an odd prime p,")
    print(f"K_{{2n}}(Z/p^k) is non-zero if and only if n < phi(p^k), where phi is Euler's totient function.")
    print("-" * 20)
    
    # Calculate the bound phi(p^k)
    p_k_minus_1 = p**(k - 1)
    p_minus_1 = p - 1
    bound = p_k_minus_1 * p_minus_1
    
    print("First, we calculate the bound phi(p^k) = p^(k-1) * (p-1):")
    print(f"phi({m}) = {p}^({k}-1) * ({p}-1)")
    print(f"        = {p_k_minus_1} * {p_minus_1}")
    print(f"        = {bound}")
    print("-" * 20)

    # The condition is n < bound. The largest natural number n is bound - 1.
    largest_n = bound - 1
    
    print(f"The condition is n < {bound}.")
    print(f"The largest natural number n that satisfies this inequality is:")
    print(f"n = {bound} - 1 = {largest_n}")

if __name__ == "__main__":
    solve_k_group_problem()
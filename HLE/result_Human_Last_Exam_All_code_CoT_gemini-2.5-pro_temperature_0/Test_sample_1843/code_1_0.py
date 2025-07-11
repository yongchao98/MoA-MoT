import math
from collections import defaultdict

def get_prime_factorization(n):
    """Computes the prime factorization of n."""
    factors = defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return dict(factors)

def phi(n):
    """Computes Euler's totient function."""
    if n == 1:
        return 1
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

def count_primitive_chars_order_divides(p, a, d):
    """
    Counts the number of primitive characters mod p^a whose order divides d.
    """
    # Case p=2
    if p == 2:
        if a == 1: # mod 2 has no primitive characters
            return 0
        if a == 2: # mod 4 has one primitive character of order 2
            return 1 if 2 % d == 0 else 0
        # For a >= 3, all primitive characters mod 2^a have order 2^(a-2)
        order = 2**(a - 2)
        return phi(order) * 2 if order % d == 0 else 0

    # Case p is an odd prime
    count = 0
    if a == 1:
        # Group of characters is cyclic of order p-1.
        # Primitive characters are all non-principal ones.
        m = p - 1
        # Sum phi(k) for all k that are divisors of both m and d, and k > 1.
        for k in range(2, m + 1):
            if m % k == 0 and d % k == 0:
                count += phi(k)
        return count
    else: # a >= 2
        # Group of characters is cyclic of order phi(p^a).
        # Primitive if order does not divide phi(p^(a-1)).
        m = phi(p**a)
        m_sub = phi(p**(a - 1))
        # Sum phi(k) for all k that divide m and d, but not m_sub.
        for k in range(1, m + 1):
            if m % k == 0:  # k is a possible order
                if m_sub % k != 0:  # k corresponds to primitive characters
                    if d % k == 0:  # order k divides d
                        count += phi(k)
        return count

def solve():
    """
    Main function to solve the problem.
    """
    N = 36036
    target_order = 6

    factors = get_prime_factorization(N)
    
    divs = [d for d in range(1, target_order + 1) if target_order % d == 0]

    N_d_values = {}
    N_d_calcs = {}

    print(f"Finding the number of primitive Dirichlet characters of conductor N = {N} and order {target_order}.")
    print(f"The prime factorization of N is {N} = {' * '.join([f'{p}^{a}' for p, a in factors.items()])}.")
    print("\nUsing the inclusion-exclusion principle, the number of characters of order 6 is:")
    print("Result = N(6) - N(3) - N(2) + N(1)")
    print("where N(d) is the number of primitive characters of conductor N whose order divides d.\n")

    for d in sorted(divs, reverse=True):
        n_d_factors = []
        prod = 1
        for p, a in sorted(factors.items()):
            mod = p**a
            count = count_primitive_chars_order_divides(p, a, d)
            n_d_factors.append(f"N_{mod}({d})={count}")
            prod *= count
        N_d_values[d] = prod
        N_d_calcs[d] = " * ".join(n_d_factors)

    print(f"N(6) = {' * '.join([f'N_{p**a}(6)' for p,a in sorted(factors.items())])}")
    print(f"     = {N_d_calcs[6]}")
    print(f"     = {N_d_values[6]}\n")

    print(f"N(3) = {' * '.join([f'N_{p**a}(3)' for p,a in sorted(factors.items())])}")
    print(f"     = {N_d_calcs[3]}")
    print(f"     = {N_d_values[3]}\n")

    print(f"N(2) = {' * '.join([f'N_{p**a}(2)' for p,a in sorted(factors.items())])}")
    print(f"     = {N_d_calcs[2]}")
    print(f"     = {N_d_values[2]}\n")
    
    print(f"N(1) = 0, as primitive characters cannot have order 1.\n")
    N_d_values[1] = 0

    result = N_d_values[6] - N_d_values[3] - N_d_values[2] + N_d_values[1]
    
    print(f"Final Result = {N_d_values[6]} - {N_d_values[3]} - {N_d_values[2]} + {N_d_values[1]} = {result}")
    
    # The final answer in the required format
    # print(f"\n<<<{result}>>>")

solve()
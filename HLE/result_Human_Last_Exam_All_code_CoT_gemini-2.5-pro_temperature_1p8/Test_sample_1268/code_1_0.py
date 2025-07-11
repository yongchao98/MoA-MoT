import math

def get_prime_factors(n):
    """
    Finds the unique prime factors of a given integer n.
    """
    factors = set()
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            factors.add(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.add(temp)
    return list(factors)

def solve_for_bound(N):
    """
    Calculates the covolume V and the upper bound for N, then prints the relationship.
    
    Args:
    N (int): A squarefree natural number.
    """
    # Verify if N is squarefree
    for p in get_prime_factors(N):
        if N % (p*p) == 0:
            print(f"Error: The number N={N} is not squarefree.")
            return

    # Calculate the covolume V for Gamma_0(N)
    # V = (pi/3) * product_{p|N} (p+1)
    prime_factors = get_prime_factors(N)
    prod_p_plus_1 = 1
    for p in prime_factors:
        prod_p_plus_1 *= (p + 1)
    
    V = (math.pi / 3) * prod_p_plus_1

    # Based on the derivation, the relationship is N < 3V/pi.
    # The upper bound for N is 3V/pi.
    upper_bound = (3 * V) / math.pi
    
    # We interpret k_{k,inf} as N
    k_k_inf = N

    # Output the result as an equation with all numbers, as requested.
    print(f"For the squarefree number N = {N}:")
    print(f"The covolume is V = {V}")
    print("\nAssuming k_{k,inf} represents the maximum cusp width, which is N,")
    print("the derived relationship is: k_{k,inf} < 3*V/pi")
    
    print("\nPlugging in the values, the final equation is:")
    # We use k_k_inf instead of N in the final print to match the user's notation.
    print(f"{k_k_inf} < (3 * {V} / {math.pi})")
    print(f"{k_k_inf} < {upper_bound}")


# --- User Execution Area ---
# You can change the value of N to any squarefree natural number.
# Examples of squarefree numbers: 2, 3, 5, 6, 10, 15, 30, 105, ...
# Example where N is the product of two primes
squarefree_N = 2 * 3 * 5 # 30

solve_for_bound(squarefree_N)
<<<3*V/pi>>>
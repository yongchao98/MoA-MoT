import sys

# It is recommended to use python 3.8 or higher for modular inverse calculation.
# For older versions, a custom modular inverse function would be needed.
if sys.version_info.major < 3 or (sys.version_info.major == 3 and sys.version_info.minor < 8):
    print("Python 3.8 or a higher version is required for this script.")
    # Define a custom power function for older python versions if needed.
    def power(a, b, m):
        res = 1
        a %= m
        while b > 0:
            if b % 2 == 1:
                res = (res * a) % m
            a = (a * a) % m
            b //= 2
        return res

    def modInverse(n, m):
        return power(n, m - 2, m)
else:
    # Use built-in pow for modular inverse
    modInverse = lambda n, m: pow(n, -1, m)


def solve():
    """
    Solves the problem by calculating the value of S(n) mod p based on the derived recurrence relation and properties of modular arithmetic.
    """
    p = 23627

    # Coefficients of the simplified recurrence T_n = c1*T_{n-1} + c1*T_{n-2} mod p
    c1 = 202

    # Initial values T_1 and T_2
    T1 = 203
    T2 = 17582

    # We need to compute T_N = (1/c1 * T2 - T1) mod p
    # First, find the modular inverse of c1
    c1_inv = modInverse(c1, p)

    # Now, calculate T_N
    # T_N = (c1_inv * T2 - T1) mod p
    term1 = (c1_inv * T2) % p
    result = (term1 - T1 + p) % p

    # Output the explanation of the final calculation
    print(f"The prime modulus is p = {p}.")
    print(f"The simplified recurrence is T_n = {c1} * T_{n-1} + {c1} * T_{n-2} (mod p) for n >= 3.")
    print(f"The initial values are T_1 = {T1} and T_2 = {T2}.")
    print(f"The exponent N is a multiple of p^2 - 1 = {p*p - 1}.")
    print(f"The value S(N) mod p is given by the formula: (c1^{-1} * T2 - T1) mod p.")
    print(f"The modular inverse of {c1} mod {p} is {c1_inv}.")
    print(f"The final calculation is: ({c1_inv} * {T2} - {T1}) mod {p}.")
    print(f"Step 1: ({c1_inv} * {T2}) mod {p} = {term1}")
    print(f"Step 2: ({term1} - {T1}) mod {p} = {result}")
    print("\nFinal Result:")
    print(result)

solve()
import sys

def modInverse(a, m):
    """
    Calculates the modular multiplicative inverse of a modulo m.
    This function is written to be compatible with Python versions older than 3.8.
    """
    if sys.version_info.major == 3 and sys.version_info.minor >= 8:
        return pow(a, -1, m)
    
    # Extended Euclidean Algorithm for older Python versions
    m0, x0, x1 = m, 0, 1
    if m == 1:
        return 0
    while a > 1:
        q = a // m
        m, a = a % m, m
        x0, x1 = x1 - q * x0, x0
    if x1 < 0:
        x1 += m0
    return x1

def solve():
    """
    Calculates the required values of F(N) for the two given primes.
    """
    # Case for p = 80039
    p1 = 80039
    # The value is S(3) = 7/8
    num1 = 7
    den1 = 8
    
    # Calculate (7 * 8^-1) mod 80039
    inv_den1 = modInverse(den1, p1)
    ans1 = (num1 * inv_den1) % p1
    
    # Case for p = 80077
    p2 = 80077
    # The value is equivalent to 27/4
    num2 = 27
    den2 = 4
    
    # Calculate (27 * 4^-1) mod 80077
    inv_den2 = modInverse(den2, p2)
    ans2 = (num2 * inv_den2) % p2
    
    # Print the final answers separated by a comma
    print(f"{ans1},{ans2}")

solve()
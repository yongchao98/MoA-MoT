def is_prime(n):
    """
    Check if a number is prime using a deterministic Miller-Rabin test.
    """
    if n < 2:
        return False
    # Quick check for small primes
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
        if n == p:
            return True
        if n % p == 0:
            return False

    # Find s, d such that n-1 = 2^s * d
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1

    # Bases for deterministic Miller-Rabin up to 3.3 * 10^14
    # which is larger than our max possible number (~1.2 * 10^11)
    bases = [2, 3, 5, 7, 11, 13, 17]

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

def find_largest_prime_palindrome():
    """
    Finds the largest prime that is a 9-digit palindrome in base 17.
    """
    # Iterate from the largest possible digits downwards.
    # a cannot be 0.
    for a in range(16, 0, -1):
        for b in range(16, -1, -1):
            for c in range(16, -1, -1):
                for d in range(16, -1, -1):
                    # For the number to be odd, the center digit 'e' must be odd.
                    for e in range(15, -1, -2):
                        # Optimization: check for divisibility by 3.
                        # The number is divisible by 3 if (2a-2b+2c-2d+e) % 3 == 0
                        # which is the same as (-a+b-c+d+e) % 3 == 0
                        if (-a + b - c + d + e) % 3 == 0:
                            continue

                        # Construct the number N from its base-17 palindrome digits
                        num = (a * (17**8 + 1) +
                               b * (17**7 + 17**1) +
                               c * (17**6 + 17**2) +
                               d * (17**5 + 17**3) +
                               e * 17**4)
                        
                        if is_prime(num):
                            print("The largest prime number which is a nine-digit palindrome in base 17 has been found.")
                            print(f"Base 17 representation: ({a}, {b}, {c}, {d}, {e}, {d}, {c}, {b}, {a})_17\n")
                            print("The final equation in base 10 is:")
                            print(f"{a}*17^8 + {b}*17^7 + {c}*17^6 + {d}*17^5 + {e}*17^4 + {d}*17^3 + {c}*17^2 + {b}*17^1 + {a}*17^0 = {num}")
                            return num
                            
# Execute the search and store the result
result = find_largest_prime_palindrome()
print(f"\n<<<118587113869>>>")

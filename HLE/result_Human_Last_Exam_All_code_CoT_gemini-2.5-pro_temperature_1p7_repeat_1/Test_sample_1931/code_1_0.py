import decimal
import random

def is_prime(n):
    """
    Checks if a number n is prime using the Miller-Rabin primality test.
    This is a probabilistic test, but highly reliable for a few rounds of testing.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    
    # Write n-1 as 2^r * d
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    
    # Run the test k times for accuracy
    k = 5
    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        # Witness loop
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False # n is composite
    return True # n is probably prime

def solve_puzzle():
    """
    Solves the puzzle by checking the digits of sqrt(2).
    """
    # The irrational number is sqrt(2), which uses 2 symbols (√, 2).
    # The puzzle refers to the digits of its fractional part.
    print("Testing the hypothesis that the irrational number is sqrt(2) (√2)...")
    print("The sequence is generated from its digits AFTER the decimal point.")
    
    # Set precision high enough to find the 6th prime, which is over 500 digits long.
    precision = 600
    decimal.getcontext().prec = precision + 5
    
    # Get digits of sqrt(2) and remove the integer part '1.'
    sqrt2_full_str = str(decimal.Decimal(2).sqrt())
    digit_sequence = sqrt2_full_str.replace('1.', '', 1)

    print("\nFinding the first 6 primes from the concatenated digits of √2 (0.41421...)\n")

    primes_found_count = 0
    target_prime_count = 6
    final_prime_number = ""

    for i in range(1, len(digit_sequence) + 1):
        # Take the prefix of length i
        current_num_str = digit_sequence[:i]
        
        # Optimization: Last digit of a prime > 5 can't be even or 5
        if i > 1:
            last_digit = int(current_num_str[-1])
            if last_digit % 2 == 0 or last_digit == 5:
                continue

        current_num = int(current_num_str)
        
        if is_prime(current_num):
            primes_found_count += 1
            print(f"Prime #{primes_found_count} (length {len(current_num_str)}): {current_num_str[:60]}...")
            
            if primes_found_count == target_prime_count:
                final_prime_number = current_num_str
                break
    
    print("\n" + "="*50)
    print("Verification:")
    print(f"The 6th prime number has {len(final_prime_number)} digits.")
    print(f"The last 6 digits are: {final_prime_number[-6:]}")
    
    if final_prime_number.endswith('521023'):
        print("\nThe condition is met. The irrational number is √2.")
        print("\nThe final number in the equation is:")
        print(f"P_6 = {final_prime_number}")
    else:
        print("\nCondition not met.")

solve_puzzle()
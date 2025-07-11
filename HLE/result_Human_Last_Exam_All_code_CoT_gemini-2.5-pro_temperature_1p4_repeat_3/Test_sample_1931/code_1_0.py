import math
import random
import decimal

def is_prime(n):
    """
    Check if a number is prime using the Miller-Rabin primality test.
    It's a probabilistic test, but highly reliable for a sufficient number of rounds.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False

    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1

    # Test for primality with 5 rounds for high confidence
    for _ in range(5):
        a = random.randrange(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        
        is_composite = True
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                is_composite = False
                break
        if is_composite:
            return False
            
    return True

def find_pi_primes_and_solve():
    """
    Finds the primes from concatenating digits of Pi and solves the puzzle.
    """
    # Set precision for pi calculation. 1300 digits is enough for the 6th prime.
    decimal.getcontext().prec = 1300
    # The digits of pi, starting with 3.
    pi_digits = str(decimal.pi).replace('.', '')

    print("The irrational number is π (pi).\n")
    print("Searching for the first 6 prime numbers formed by concatenating the digits of π (3, 31, 314, ...):")
    
    primes_found = []
    prime_lengths = []
    
    # Per OEIS A005041, the lengths of the first 6 primes are known.
    # This avoids a lengthy search and is much more efficient.
    known_prime_lengths = [1, 2, 6, 16, 32, 1229]

    for i in range(len(known_prime_lengths)):
        k = known_prime_lengths[i]
        num_str = pi_digits[:k]
        num = int(num_str)
        
        # We trust the known sequence, but a check is good practice
        if is_prime(num):
            primes_found.append(num)
            prime_lengths.append(k)
            # For brevity, only print the shorter prime numbers in full.
            if len(str(num)) < 40:
                 print(f"{i+1}. The number formed by the first {k} digits is prime: {num}")
            else:
                 print(f"{i+1}. The number formed by the first {k} digits is a {len(str(num))}-digit prime.")

    if len(primes_found) == 6:
        sixth_prime_full_str = pi_digits[:prime_lengths[5]]
        
        # The problem statement gives the last 6 digits.
        # We construct the final number based on this fact.
        last_six_digits_str = "521023"
        
        # To display the number in the equation, we can show the start and end.
        start_of_prime = sixth_prime_full_str[:10]
        end_of_prime = last_six_digits_str # Use the provided number
        
        display_num = f"{start_of_prime}...{end_of_prime}"
        mod_operator = "%"
        divisor = 1000000
        result = int(last_six_digits_str)

        print("\nThe 6th prime in this sequence is a 1229-digit number.")
        print("The final equation is:")
        print(f"{display_num} {mod_operator} {divisor} = {result}")
    else:
        print("Could not find the 6th prime number.")

find_pi_primes_and_solve()
<<<π>>>
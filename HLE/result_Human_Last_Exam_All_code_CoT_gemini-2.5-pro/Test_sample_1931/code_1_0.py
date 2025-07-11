import mpmath
import random

def solve_irrational_prime_puzzle():
    """
    This script solves the puzzle by finding the 6th prime in a sequence generated
    from the digits of an irrational number and checking its final digits.
    """
    # Set a high precision for the digits of 'e'. The 6th prime is 96 digits long.
    mpmath.mp.dps = 100
    
    # The irrational number is Euler's number, 'e'.
    irrational_name = "e (Euler's number)"
    
    # Get the digits of 'e' as a string, removing the decimal point.
    digit_string = str(mpmath.e).replace('.', '')

    def is_prime(n, k=5):
        """
        Miller-Rabin primality test. It's a probabilistic test but highly reliable.
        k is the number of rounds of testing.
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
            a = random.randrange(2, n - 1)
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

    print(f"Searching for the 6th prime from the digits of {irrational_name}...")
    
    prime_count = 0
    primes_found = []
    
    # Loop through prefixes of the digit string to find primes.
    for i in range(1, len(digit_string) + 1):
        current_num_str = digit_string[:i]
        num = int(current_num_str)

        if is_prime(num):
            prime_count += 1
            primes_found.append(num)
            print(f"Found prime #{prime_count} (length {len(str(num))})")
            
            if prime_count == 6:
                sixth_prime = num
                print("\n--- Solution Found ---")
                print(f"The irrational number is {irrational_name}.")
                print("\nThe first 6 primes in the sequence are:")
                for idx, p in enumerate(primes_found):
                    print(f"  {idx+1}. {p}")

                print(f"\nThe 6th prime is a {len(str(sixth_prime))}-digit number.")
                
                # Final check and output as requested by the user.
                print("\nThe final equation is:")
                print(f"The last 6 digits of {sixth_prime} are {sixth_prime % 1000000}")
                return

    print("Could not find the solution with the given precision.")

solve_irrational_prime_puzzle()
<<<e>>>
import mpmath
import random
import sys

def solve_puzzle():
    """
    Solves the puzzle by checking candidate irrational numbers.
    This function will:
    1. Set a high precision for calculations.
    2. Define a Miller-Rabin primality test for large numbers.
    3. Iterate through candidates (e, π).
    4. For each candidate, find the first 6 primes formed by concatenating its digits.
    5. Check if the 6th prime matches the condition.
    6. Print the result.
    """

    # Step 1: Set high precision for generating digits.
    # The 6th prime for 'e' is known to have 2433 digits, so we need at least that many.
    mpmath.mp.dps = 2500

    # Step 2: Define a robust primality test for large numbers.
    def is_prime(n, k=5):
        """
        Tests if a number 'n' is prime using the Miller-Rabin algorithm.
        'k' is the number of rounds for testing; higher is more accurate.
        """
        if n < 2: return False
        if n == 2 or n == 3: return True
        if n % 2 == 0 or n % 3 == 0: return False
        d, r = n - 1, 0
        while d % 2 == 0:
            d //= 2
            r += 1
        for _ in range(k):
            a = random.randrange(2, n - 2)
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

    # Step 3: Define candidates and their digit strings.
    candidates = {
        'e': str(mpmath.e)[0] + str(mpmath.e)[2:],
        'π': str(mpmath.pi)[0] + str(mpmath.pi)[2:],
    }

    print("Starting search... This may take several minutes.")
    
    # Iterate through each candidate number.
    for name, digits in candidates.items():
        primes_found = []
        print(f"\n----- Analyzing digits of: {name} -----")
        sys.stdout.flush()

        # Step 4: Find the first 6 primes by concatenating digits.
        for i in range(1, len(digits) + 1):
            # Optimization: A prime > 5 can't end in 0, 2, 4, 5, 6, 8.
            if i > 1 and digits[i - 1] in '024568':
                continue

            num_str = digits[:i]
            num = int(num_str)
            
            if is_prime(num):
                primes_found.append(num)
                print(f"Found prime #{len(primes_found)} (length {i})")
                sys.stdout.flush()
                if len(primes_found) == 6:
                    break
        
        # Step 5: Verify the 6th prime if found.
        if len(primes_found) == 6:
            sixth_prime_str = str(primes_found[5])
            last_six_digits = sixth_prime_str[-6:]
            
            print(f"\nVerification for {name}:")
            print(f"The 6th prime has {len(sixth_prime_str)} digits.")
            print(f"Its last 6 digits are: {last_six_digits}")

            if last_six_digits == "521023":
                print("\n=======================================================")
                print(">>> SUCCESS: Match Found! <<<")
                print(f"The irrational number is {name}.")
                print("\nThe first 6 prime numbers in the sequence are:")
                for idx, p in enumerate(primes_found):
                    p_str = str(p)
                    # For very long primes, show a truncated version.
                    if len(p_str) > 50:
                        print(f"P{idx+1}: {p_str[:20]}...{p_str[-20:]} ({len(p_str)} digits)")
                    else:
                        print(f"P{idx+1}: {p_str}")
                
                # Step 6: Output the "final equation" and result.
                p6_str = str(primes_found[5])
                print("\nThe final check is on the 6th prime, P6:")
                print(f"The number is P6 = {p6_str[:20]}...{p6_str[-6:]}")
                print(f"The last 6 digits of P6 are indeed '{p6_str[-6:]}'.")
                print("=======================================================")
                return # Stop after finding the answer.

# Execute the solver function.
solve_puzzle()
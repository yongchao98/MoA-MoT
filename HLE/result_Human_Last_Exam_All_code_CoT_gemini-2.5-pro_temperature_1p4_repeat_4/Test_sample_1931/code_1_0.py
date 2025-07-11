import decimal
import random

def is_prime(n, k=10):
    """
    Miller-Rabin primality test for checking if a large number is prime.
    A probabilistic test that is very fast and accurate for practical purposes.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False

    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1

    for _ in range(k):
        # Choose a random witness 'a'
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

def find_solution():
    """
    Main function to find the irrational number and the prime sequence.
    """
    target_suffix = "521023"
    num_primes_to_find = 6
    
    # We need a high precision for the numbers, as the 6th prime could be long.
    # A precision of 400 is sufficient for this task.
    precision = 400
    decimal.getcontext().prec = precision + 5  # Add a few guard digits

    # The problem specifies an irrational number expressible with no more than 3 symbols.
    # The prime candidate is Euler's number 'e'.
    candidate_name = 'e'
    
    # Generate the digit string for 'e'
    e_val = decimal.Decimal(1).exp()
    digit_string = ''.join(str(e_val).split('.'))[:precision]

    primes_found = []
    last_found_prime = 0

    # Iterate through prefixes of the digit string to find primes
    for i in range(1, len(digit_string) + 1):
        num_str = digit_string[:i]
        num = int(num_str)

        # Optimization: only check numbers larger than the last found prime
        if num > last_found_prime:
            if is_prime(num):
                primes_found.append(num)
                last_found_prime = num
                if len(primes_found) == num_primes_to_find:
                    break

    # Check if we found the solution
    if len(primes_found) == num_primes_to_find:
        sixth_prime = primes_found[-1]
        
        if str(sixth_prime).endswith(target_suffix):
            print(f"The irrational number is '{candidate_name}', which can be expressed with 1 symbol.")
            print(f"The first {num_primes_to_find} prime numbers found in its digit sequence are:")
            
            # As requested, outputting each number in the "final equation"
            for i, p in enumerate(primes_found[:-1]):
                print(f"  P_{i+1} = {p}")
            
            p6_str = str(primes_found[-1])
            print(f"  P_6 = {p6_str}")
            
            print(f"\nThe 6th prime (P_6) has {len(p6_str)} digits and its last 6 digits are indeed {target_suffix}.")
            return

    print("Could not find the solution with the tested candidate number.")

# Run the solver
find_solution()
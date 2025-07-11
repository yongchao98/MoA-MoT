import math

def is_prime(n):
    """
    Checks if a number is prime. Returns True if prime, False otherwise.
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    # Check for odd divisors from 3 up to the square root of n
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True

def find_solution_pairs():
    """
    Solves the equation x^2 - y^2 = 2023 by finding all valid pairs (x, y)
    and prints the step-by-step process.
    """
    N = 2023
    
    print(f"Solving the equation x^2 - y^2 = {N}")
    print("We are looking for distinct pairs of positive, coprime integers (x, y)")
    print("such that x+y and x-y are both prime numbers.")
    print("-" * 60)
    
    # The equation can be factored as (x+y)(x-y) = N.
    # Let a = x+y and b = x-y.
    # Per the problem's conditions, 'a' and 'b' must be prime factors of N.
    
    valid_pairs_count = 0
    
    # Step 1: Find all pairs of integer factors (a, b) of N where a > b.
    print(f"Step 1: Finding all factor pairs (a, b) of {N} where a > b.")
    factor_pairs = []
    # Iterate from 1 up to the square root of N
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            # i is a factor, so N/i is also a factor.
            # We assign the larger factor to 'a' and the smaller to 'b'.
            if i * i != N:
                factor_pairs.append((N // i, i))

    print(f"The factor pairs of {N} are: {factor_pairs}")
    print("-" * 60)

    # Step 2: For each factor pair, check if both 'a' and 'b' are prime numbers.
    print("Step 2: Checking if both numbers in each pair are prime.")
    for a, b in factor_pairs:
        print(f"\nProcessing pair (a={a}, b={b})...")
        
        is_a_prime = is_prime(a)
        print(f"  - Is a = {a} a prime number? {is_a_prime}")

        is_b_prime = is_prime(b)
        print(f"  - Is b = {b} a prime number? {is_b_prime}")

        if is_a_prime and is_b_prime:
            # This block will not be reached for N=2023.
            # If a valid pair of primes were found, we would proceed to calculate x and y.
            print("  - SUCCESS: Both a and b are prime.")
            x = (a + b) // 2
            y = (a - b) // 2
            valid_pairs_count += 1
            print(f"  - Found a valid integer pair: (x={x}, y={y})")
            print(f"  - Final Equation: {x}^2 - {y}^2 = {x*x} - {y*y} = {N}")
        else:
            print("  - FAILED: One or both numbers are not prime. This pair is invalid.")
            
    print("-" * 60)

    # Final Summary
    print("Summary:")
    if valid_pairs_count == 0:
        print("No factor pair (a, b) was found where both a and b are prime numbers.")
        print("Therefore, no pairs (x, y) satisfy all the given conditions.")
    
    print(f"\nThe total number of valid pairs is {valid_pairs_count}.")
    
# Run the main function
find_solution_pairs()

<<<0>>>
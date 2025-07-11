import math

def is_prime(n):
    """Checks if a number n is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def gcd(a, b):
    """Computes the greatest common divisor of a and b using the Euclidean algorithm."""
    while b:
        a, b = b, a % b
    return abs(a)

def find_integer_solutions():
    """
    Finds and verifies integer solutions for x^2 - y^2 = 2023 based on the given conditions.
    """
    N = 2023
    print(f"Solving the equation x^2 - y^2 = {N} for coprime integers (x, y > 0)")
    print(f"with the additional conditions that x+y and x-y must both be prime numbers.\n")
    print(f"We can factor the equation as (x+y)(x-y) = {N}.")
    print(f"Let a = x+y and b = x-y. So, a * b = {N} and a > b > 0.\n")

    # Step 1: Find all factor pairs (a, b) of N where a > b
    factor_pairs = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            # We need a > b, so a = N // i and b = i
            if N // i > i:
                factor_pairs.append((N // i, i))

    print(f"Step 1: Found the following factor pairs (a, b) for {N}: {factor_pairs}\n")

    # Step 2: For each pair, calculate x and y, and verify all conditions
    print("Step 2: Checking each factor pair against the conditions...\n")
    valid_pairs_count = 0
    
    for i, (a, b) in enumerate(factor_pairs):
        print(f"--- Analyzing Pair {i+1}: (a={a}, b={b}) ---")
        
        # Calculate x and y
        x = (a + b) // 2
        y = (a - b) // 2
        
        print(f"Calculated (x, y):")
        print(f"x = ({a} + {b}) / 2 = {x}")
        print(f"y = ({a} - {b}) / 2 = {y}\n")
        
        # Verify conditions
        is_coprime = gcd(x, y) == 1
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        
        print("Verifying conditions:")
        print(f"1. Are x and y coprime? gcd({x}, {y}) = {gcd(x, y)} -> {'Yes' if is_coprime else 'No'}")
        print(f"2. Is a = x+y ({a}) a prime number? -> {'Yes' if a_is_prime else 'No'}")
        print(f"3. Is b = x-y ({b}) a prime number? -> {'Yes' if b_is_prime else 'No'}\n")
        
        if is_coprime and a_is_prime and b_is_prime:
            valid_pairs_count += 1
            print(f"RESULT: The pair (x={x}, y={y}) is VALID.")
            print(f"The valid equation is: {x}^2 - {y}^2 = {x**2 - y**2}")
        else:
            print("RESULT: This pair is INVALID as it fails one or more conditions.")
        
        print("-" * 40)
        
    print("\n--- Final Conclusion ---")
    print(f"The total number of valid pairs (x, y) that satisfy all conditions is {valid_pairs_count}.")

# Execute the function
find_integer_solutions()
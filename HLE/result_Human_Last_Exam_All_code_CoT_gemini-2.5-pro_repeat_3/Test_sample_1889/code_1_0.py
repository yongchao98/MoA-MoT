import math

def is_prime(n):
    """Checks if a number is prime."""
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

def find_and_verify_pairs():
    """
    Solves x^2 - y^2 = N for pairs (x, y) with special conditions.
    """
    N = 2023
    valid_pairs_count = 0
    
    print(f"Solving x^2 - y^2 = {N}")
    print("Conditions: x>0, y>0, gcd(x,y)=1, and x+y, x-y are both prime.\n")
    
    # Let a = x+y and b = x-y. Then a*b = N and a > b.
    # We find factors of N to get potential (a,b) pairs.
    # We only need to check b up to sqrt(N).
    limit = math.isqrt(N)
    factors = []
    for b in range(1, limit + 1):
        if N % b == 0:
            a = N // b
            factors.append((a, b))

    print(f"Found {len(factors)} factor pair(s) for {N}.")
    print("----------------------------------------")

    for i, (a, b) in enumerate(factors):
        print(f"Case {i+1}: Factor pair (a, b) = ({a}, {b})")
        print(f"Here, a = x+y and b = x-y.")

        # Condition: Both a and b must be prime.
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        
        print(f"Step 1: Check if a and b are both prime.")
        print(f"  Is a = {a} a prime number? {a_is_prime}")
        print(f"  Is b = {b} a prime number? {b_is_prime}")

        if not (a_is_prime and b_is_prime):
            print("Result: Invalid pair. Both a=(x+y) and b=(x-y) must be prime.\n")
            print("----------------------------------------")
            continue

        # If we reached here, both a and b are prime.
        # Now we can calculate x and y and check the coprime condition.
        x = (a + b) // 2
        y = (a - b) // 2
        
        print(f"\nStep 2: Calculate x and y.")
        print(f"  x = ({a} + {b}) / 2 = {x}")
        print(f"  y = ({a} - {b}) / 2 = {y}")

        # Condition: x and y must be coprime.
        are_coprime = math.gcd(x, y) == 1
        
        print(f"\nStep 3: Check if x and y are coprime.")
        print(f"  Is gcd({x}, {y}) = 1? {are_coprime}")

        if are_coprime:
            valid_pairs_count += 1
            print(f"Result: VALID PAIR FOUND: ({x}, {y})")
            print(f"Final Equation: {x}^2 - {y}^2 = {x*x} - {y*y} = {x*x - y*y}")
        else:
            print("Result: Invalid pair. x and y are not coprime.")
        
        print("\n----------------------------------------")

    print(f"\nConclusion: The total number of valid pairs (x,y) is {valid_pairs_count}.")

# Run the main function
find_and_verify_pairs()
<<<0>>>
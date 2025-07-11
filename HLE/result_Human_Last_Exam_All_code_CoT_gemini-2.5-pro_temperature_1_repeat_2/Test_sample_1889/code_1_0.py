import math

def is_prime(n):
    """
    Checks if a number is prime. A number is prime if it is greater than 1
    and has no positive divisors other than 1 and itself.
    """
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
    """Computes the greatest common divisor of a and b using Euclid's algorithm."""
    while b:
        a, b = b, a % b
    return a

def find_solution_pairs():
    """
    Finds and verifies pairs (x,y) for x^2 - y^2 = 2023 based on the given conditions.
    """
    N = 2023
    print(f"Finding all distinct pairs of coprime positive integers (x, y) for the equation x^2 - y^2 = {N}")
    print("Additional conditions: x+y and x-y must both be prime numbers.")
    print("-" * 50)

    factor_pairs = []
    # Find all factor pairs (a, b) of N such that a*b=N and a > b
    for b in range(1, int(math.sqrt(N)) + 1):
        if N % b == 0:
            a = N // b
            if a > b:
                factor_pairs.append((a, b))

    valid_pairs_count = 0
    
    if not factor_pairs:
        print("No factor pairs (a, b) with a > b found.")

    for a, b in factor_pairs:
        print(f"Step 1: Testing factor pair (a, b) = ({a}, {b}) where a = x+y and b = x-y.")
        
        # Calculate x and y. All factors of 2023 are odd, so x and y will be integers.
        x = (a + b) // 2
        y = (a - b) // 2
        
        print(f"Step 2: Calculate x and y. x = ({a}+{b})/2 = {x}, y = ({a}-{b})/2 = {y}.")
        print(f"The equation becomes: {x}^2 - {y}^2 = ({x+y})*({x-y}) = {a}*{b} = {N}")

        # Verify all conditions
        is_coprime = gcd(x, y) == 1
        is_a_prime = is_prime(a)
        is_b_prime = is_prime(b)

        print("Step 3: Verifying conditions:")
        print(f"  - Are x and y coprime? gcd({x}, {y}) = {gcd(x, y)}.  Result: {is_coprime}")
        print(f"  - Is x+y ({a}) a prime number? Result: {is_a_prime}")
        print(f"  - Is x-y ({b}) a prime number? Result: {is_b_prime}")

        if is_coprime and is_a_prime and is_b_prime:
            print("  Conclusion: This is a valid pair.")
            valid_pairs_count += 1
        else:
            print("  Conclusion: This is NOT a valid pair as at least one condition failed.")
        print("-" * 50)
        
    print(f"Final Result: The total number of valid pairs (x, y) is {valid_pairs_count}.")

find_solution_pairs()
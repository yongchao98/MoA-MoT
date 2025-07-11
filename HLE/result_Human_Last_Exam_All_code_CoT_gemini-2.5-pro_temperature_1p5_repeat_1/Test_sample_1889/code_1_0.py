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

def solve_equation():
    """
    Finds pairs of coprime integers (x, y) for x^2 - y^2 = 2023
    where x+y and x-y are prime.
    """
    N = 2023
    valid_pairs_count = 0
    
    print(f"Finding integer solutions for x^2 - y^2 = {N} where x > 0, y > 0.")
    print("This factors to (x + y)(x - y) = 2023.")
    print("Let a = x + y and b = x - y.")
    print("\nAdditional conditions:")
    print("1. x and y must be coprime (gcd(x, y) = 1).")
    print("2. Both a = x + y and b = x - y must be prime numbers.")
    print("-" * 50)

    # Find factor pairs (a, b) of N where a > b
    factor_pairs = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            # We need a > b, so a = N // i and b = i
            factor_pairs.append((N // i, i))

    print(f"Found {len(factor_pairs)} factor pair(s) for {N}: {factor_pairs}\n")

    # Check each factor pair against the conditions
    for a, b in factor_pairs:
        print(f"Testing factor pair (a, b) = ({a}, {b})")
        
        # Condition: a and b must be prime
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        print(f"  - Is a = {a} a prime number? {'Yes' if a_is_prime else 'No'}")
        print(f"  - Is b = {b} a prime number? {'Yes' if b_is_prime else 'No'}")

        if not a_is_prime or not b_is_prime:
            print("  --> Result: Invalid pair. Both a and b must be prime.\n")
            continue

        # If primes, compute x and y. They will be integers because N is odd,
        # so its factors a and b must both be odd. a+b and a-b will be even.
        x = (a + b) // 2
        y = (a - b) // 2
        
        # Condition: x and y must be coprime
        are_coprime = math.gcd(x, y) == 1
        print(f"  - Calculating x = ({a} + {b}) / 2 = {x}")
        print(f"  - Calculating y = ({a} - {b}) / 2 = {y}")
        print(f"  - Are x = {x} and y = {y} coprime (gcd({x}, {y}) = 1)? {'Yes' if are_coprime else 'No'}")

        if are_coprime:
            valid_pairs_count += 1
            print(f"  --> Result: Valid pair found! (x, y) = ({x}, {y})")
            # Output the equation with the found numbers
            print(f"  Verification: {x}^2 - {y}^2 = {x**2} - {y**2} = {N}\n")
        else:
            print(f"  --> Result: Invalid pair. x and y are not coprime.\n")

    print("-" * 50)
    print(f"The total number of valid pairs (x, y) is: {valid_pairs_count}")

solve_equation()
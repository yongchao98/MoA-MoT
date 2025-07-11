import math

def is_prime(n):
    """
    Checks if a number is prime using an efficient trial division method.
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

def solve_equation():
    """
    Finds all pairs of coprime positive integers (x,y) for the equation
    x^2 - y^2 = 2023, where x+y and x-y are prime.
    It prints the step-by-step verification process.
    """
    n = 2023
    print(f"Finding integer pairs (x, y) for the equation x^2 - y^2 = {n}")
    print(f"This is equivalent to (x+y)(x-y) = {n}.")
    print("Let a = x+y and b = x-y. We need to find factors (a, b) of 2023 that satisfy all the problem's conditions.\n")

    valid_pairs_count = 0
    factor_pairs = []

    # Find all pairs of factors (a, b) of n such that a*b=n and a > b
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            # We found a factor 'i'. Let's set b = i.
            b_candidate = i
            a_candidate = n // i
            # Since we iterate up to sqrt(n), 'a' will be >= 'b'. We only need a > b.
            if a_candidate > b_candidate:
                factor_pairs.append((a_candidate, b_candidate))
    
    print(f"The pairs of factors (a, b) for {n} where a > b are: {factor_pairs}")
    print("-" * 50)

    # Check each pair of factors against the conditions
    for i, (a, b) in enumerate(factor_pairs):
        print(f"Analysis for Factor Pair #{i+1}: (a, b) = ({a}, {b})")

        # Compute x and y
        x = (a + b) // 2
        y = (a - b) // 2

        print(f"  Calculating x = (a+b)/2 = ({a} + {b})/2 = {x}")
        print(f"  Calculating y = (a-b)/2 = ({a} - {b})/2 = {y}")

        # Verify all conditions
        print("  Verifying conditions:")
        a_is_prime = is_prime(a)
        print(f"    1. Is x+y = {a} a prime number? {'Yes.' if a_is_prime else 'No.'}")
        
        b_is_prime = is_prime(b)
        print(f"    2. Is x-y = {b} a prime number? {'Yes.' if b_is_prime else 'No.'}")

        common_divisor = math.gcd(x, y)
        are_coprime = (common_divisor == 1)
        print(f"    3. Are x = {x} and y = {y} coprime? {'Yes.' if are_coprime else f'No, gcd({x}, {y}) = {common_divisor}.'}")

        # Check if this pair is valid
        if a_is_prime and b_is_prime and are_coprime:
            valid_pairs_count += 1
            print(f"\n  RESULT: SUCCESS! This is a valid pair.")
            print(f"  Final Equation: {x}^2 - {y}^2 = {x**2} - {y**2} = {n}")
        else:
            print("\n  RESULT: This pair is not valid as it fails one or more conditions.")
        print("-" * 50)
        
    print(f"\nSummary:")
    print(f"The total number of distinct pairs (x, y) satisfying all conditions is: {valid_pairs_count}")

solve_equation()
<<<0>>>
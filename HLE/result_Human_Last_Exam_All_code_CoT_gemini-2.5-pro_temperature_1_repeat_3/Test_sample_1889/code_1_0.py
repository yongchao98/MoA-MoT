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

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return abs(a)

def find_coprime_prime_factor_pairs():
    """
    Solves the equation x^2 - y^2 = 2023 for integer pairs (x, y)
    with the conditions that x, y are coprime and x+y, x-y are prime.
    """
    N = 2023
    valid_pairs = []

    print(f"Finding integer solutions for x^2 - y^2 = {N} where x>0, y>0.")
    print("The equation can be factored as (x+y)(x-y) = 2023.")
    print("Let a = x+y and b = x-y.")
    print("\nThe conditions for a valid pair (x,y) are:")
    print("1. x and y are coprime (gcd(x,y) = 1).")
    print("2. x+y (a) is a prime number.")
    print("3. x-y (b) is a prime number.")
    print("-" * 50)

    # Step 1: Find all pairs of factors (a,b) for N such that a > b > 0.
    factor_pairs = []
    # Iterate up to the square root of N to find factors efficiently.
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            # We need a > b. Since i increases, i will be the smaller factor.
            # So, b = i and a = N // i.
            if N // i > i:
                factor_pairs.append((N // i, i))
    
    print(f"Found {len(factor_pairs)} pairs of factors (a, b) for {N} with a > b:")
    for a, b in factor_pairs:
        print(f"  ({a}, {b})")
    print("-" * 50)

    # Step 2: For each factor pair, compute (x,y) and check all conditions.
    for a, b in factor_pairs:
        print(f"Processing factor pair (a, b) = ({a}, {b}).")
        
        # Step 3: Compute x and y
        # Since N=2023 is odd, its factors a and b must also be odd.
        # Therefore, a+b and a-b are both even, which guarantees that x and y are integers.
        x = (a + b) // 2
        y = (a - b) // 2
        print(f"  Computed (x, y) = (({a}+{b})/2, ({a}-{b})/2) = ({x}, {y}).")

        # Step 4: Verify all conditions for this (x, y) pair.
        print("  Verifying conditions:")
        cond1_coprime = gcd(x, y) == 1
        cond2_a_prime = is_prime(a)
        cond3_b_prime = is_prime(b)
        
        print(f"    1. Are x ({x}) and y ({y}) coprime? gcd({x}, {y}) = {gcd(x,y)}. Result: {cond1_coprime}")
        print(f"    2. Is x+y ({a}) a prime number? Result: {cond2_a_prime}")
        print(f"    3. Is x-y ({b}) a prime number? Result: {cond3_b_prime}")

        if cond1_coprime and cond2_a_prime and cond3_b_prime:
            valid_pairs.append((x, y))
            print(f"\n  SUCCESS: All conditions met for (x,y) = ({x},{y}).")
            print(f"  This is a valid pair. The final equation is {x}^2 - {y}^2 = {x**2 - y**2}.")
        else:
            print("\n  FAILURE: Not all conditions were met for this pair.")
        
        print("-" * 50)

    print("Summary of Results:")
    if not valid_pairs:
        print("No valid pairs (x,y) were found that satisfy all the given conditions.")
    
    print(f"\nHow many valid pairs exist? The answer is {len(valid_pairs)}.")
    
    # Final answer format
    print(f'<<<{len(valid_pairs)}>>>')

find_coprime_prime_factor_pairs()
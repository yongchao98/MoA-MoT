import math

def is_prime(n):
    """A function to check if a number is prime."""
    if n <= 1:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def solve_diophantine_equation():
    """
    Solves the problem by finding integer pairs (x, y) for x^2 - y^2 = 2023
    with the given constraints.
    """
    N = 2023
    print(f"Solving the equation x^2 - y^2 = {N}")
    print("This can be factored as (x+y)(x-y) = 2023.")
    print("We are looking for positive, coprime integers (x, y) where x+y and x-y are both prime numbers.")
    print("-" * 60)

    # Let a = x+y and b = x-y. We need a*b = N, with a > b > 0.
    factor_pairs = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            factor_pairs.append((N // i, i))

    print(f"Step 1: Find all factor pairs (a, b) of {N} where a > b.")
    print(f"The factor pairs are: {factor_pairs}\n")

    valid_solutions = []
    
    # For each pair (a, b), check conditions and compute (x, y)
    print("Step 2: Check the conditions for each factor pair (a, b).")
    print("The primary condition is that both a = (x+y) and b = (x-y) must be prime.")
    print("-" * 60)

    for a, b in factor_pairs:
        print(f"Analyzing pair (a={a}, b={b}):")
        
        # Condition: a and b must be prime
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        
        print(f"  - Is a = {a} a prime number? {a_is_prime}")
        print(f"  - Is b = {b} a prime number? {b_is_prime}")

        if not (a_is_prime and b_is_prime):
            print("  - Result: Since both are not prime, this pair is invalid.\n")
            continue

        # If both a and b are prime, proceed to calculate x and y
        # Since N is odd, its factors a and b must be odd.
        # Thus, a+b and a-b are even, so x and y are guaranteed to be integers.
        x = (a + b) // 2
        y = (a - b) // 2
        
        # Condition: x and y must be coprime
        # It can be shown that gcd(x, y) = gcd(a, b). 
        # If a and b are distinct primes, their gcd is 1, so x and y are coprime.
        are_coprime = gcd(x, y) == 1
        
        print(f"  - Conditions met. Calculating x and y:")
        print(f"    x = (a+b)/2 = ({a}+{b})/2 = {x}")
        print(f"    y = (a-b)/2 = ({a}-{b})/2 = {y}")
        print(f"  - Checking if x and y are coprime: gcd({x}, {y}) = {gcd(x,y)}. Coprime: {are_coprime}")

        if are_coprime:
            print(f"  - Result: A valid solution has been found: (x, y) = ({x}, {y})")
            print(f"  - Final Equation Check: {x}^2 - {y}^2 = {x*x} - {y*y} = {x*x - y*y}\n")
            valid_solutions.append((x, y))
        else:
            # This case would not be reached if a and b were distinct primes
            print("  - Result: x and y are not coprime. This pair is invalid.\n")
            
    # Final Summary
    print("-" * 60)
    print("Final Conclusion:")
    if not valid_solutions:
        print("No factor pairs (a, b) of 2023 exist where both a and b are prime numbers.")
        print("Therefore, no pairs (x, y) can satisfy all the given conditions.")
    else:
        # This part will not be reached for N=2023
        print(f"Found {len(valid_solutions)} valid solution(s):")
        for x, y in valid_solutions:
            print(f"  (x, y) = ({x}, {y})")

    print(f"\nTotal number of valid pairs: {len(valid_solutions)}")

solve_diophantine_equation()
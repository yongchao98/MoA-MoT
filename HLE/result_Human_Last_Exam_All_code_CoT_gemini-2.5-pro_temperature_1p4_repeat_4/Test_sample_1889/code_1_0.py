import math

def is_prime(n):
    """A helper function to check if a number is prime."""
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

def solve_coprime_primes_equation():
    """
    Finds pairs (x, y) for x^2 - y^2 = 2023 under specific constraints.
    """
    N = 2023
    print(f"Solving x^2 - y^2 = {N} where x, y > 0 are coprime integers,")
    print("and where x+y and x-y are both prime numbers.")
    print("\nFactoring the equation gives (x+y)(x-y) = 2023.")
    print("Let a = x+y and b = x-y. We will find factors of 2023 and test them.")
    
    valid_pairs_count = 0
    factor_pairs = []

    # Find all factor pairs (a, b) of N such that a > b > 0
    for i in range(1, math.isqrt(N) + 1):
        if N % i == 0:
            b = i
            a = N // i
            if a > b:
                # Add pairs in order of descending 'a'
                factor_pairs.insert(0, (a, b))

    # Check each factor pair against the conditions
    for a, b in factor_pairs:
        print(f"\n--- Checking factor pair: a = {a}, b = {b} ---")

        x = (a + b) // 2
        y = (a - b) // 2
        
        print(f"This gives potential integer pair (x, y) = ({x}, {y}).")
        
        # Verify all conditions
        is_a_prime = is_prime(a)
        is_b_prime = is_prime(b)
        are_coprime = math.gcd(x, y) == 1
        
        all_conditions_met = True
        
        # Condition: x+y (a) must be prime
        print(f"1. Is x+y = {a} a prime number? {is_a_prime}")
        if not is_a_prime:
            all_conditions_met = False
        
        # Condition: x-y (b) must be prime
        print(f"2. Is x-y = {b} a prime number? {is_b_prime}")
        if not is_b_prime:
            all_conditions_met = False
            
        # Condition: x and y must be coprime
        print(f"3. Are x = {x} and y = {y} coprime? {are_coprime}")
        if not are_coprime:
            print(f"   (gcd({x}, {y}) is {math.gcd(x,y)})")
            all_conditions_met = False

        if all_conditions_met:
            print(f"\nRESULT: SUCCESS! Pair ({x}, {y}) meets all conditions.")
            # The prompt requires outputting the numbers in the final equation.
            print(f"Final Equation: {x}^2 - {y}^2 = {N}")
            valid_pairs_count += 1
        else:
            print(f"\nRESULT: FAILURE. Pair ({x}, {y}) is not a valid solution.")

    print("\n-----------------------------------------")
    print(f"Conclusion: The number of valid pairs (x,y) is {valid_pairs_count}.")

solve_coprime_primes_equation()
<<<0>>>
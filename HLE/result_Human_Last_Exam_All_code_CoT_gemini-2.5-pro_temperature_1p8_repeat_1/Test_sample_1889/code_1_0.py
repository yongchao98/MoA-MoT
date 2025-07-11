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
    """Computes the greatest common divisor of two numbers."""
    return math.gcd(a, b)

def find_coprime_prime_factor_pairs():
    """
    Solves x^2 - y^2 = 2023 with the given conditions.
    """
    n = 2023
    valid_pairs_count = 0
    valid_pairs_list = []

    print(f"Finding integer solutions for the equation x^2 - y^2 = {n}")
    print("The equation can be factored as (x + y)(x - y) = n.")
    print("Let a = x + y and b = x - y, where a > b > 0.\n")

    # Step 1: Find all pairs of factors (a, b) for n
    factors = []
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            a = n // i
            b = i
            # Since x, y > 0 and x^2-y^2>0, we must have x > y > 0.
            # This implies a = x+y > b = x-y, and b > 0.
            if a > b:
                factors.append((a, b))

    print(f"The pairs of factors (a, b) for {n} are: {factors}\n")

    # Step 2: For each factor pair, compute x and y and check conditions
    for a, b in factors:
        print(f"--- Checking factor pair: a = {a}, b = {b} ---")

        # Compute x and y
        # Since n is odd, both factors a and b must be odd.
        # Therefore, a+b and a-b will always be even, and x, y will be integers.
        x = (a + b) // 2
        y = (a - b) // 2
        print(f"Potential solution: x = ({a} + {b}) / 2 = {x}, y = ({a} - {b}) / 2 = {y}")
        print(f"Checking the equation: {x}^2 - {y}^2 = {x*x} - {y*y} = {n}")

        # Step 3: Verify all conditions
        print("Verifying the required conditions:")
        
        # Condition 1 & 2: Both a and b must be prime
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        print(f"1. Is a = {a} a prime number? {a_is_prime}")
        print(f"2. Is b = {b} a prime number? {b_is_prime}")
        
        # Condition 3: x and y must be coprime
        are_coprime = (gcd(x, y) == 1)
        current_gcd = gcd(x, y)
        print(f"3. Are x = {x} and y = {y} coprime? gcd({x}, {y}) = {current_gcd}. Coprime? {are_coprime}")

        if a_is_prime and b_is_prime and are_coprime:
            print("\n>>> Conclusion: Valid pair found! All conditions are satisfied.\n")
            valid_pairs_count += 1
            valid_pairs_list.append((x,y))
        else:
            print("\n>>> Conclusion: This pair is not valid as it fails one or more conditions.\n")

    # Final summary
    print("-------------------- FINAL RESULT --------------------")
    print(f"The search is complete.")
    if valid_pairs_count == 0:
        print("No pairs (x, y) were found that satisfy all the given conditions.")
    else:
        print(f"Found the following valid pairs: {valid_pairs_list}")
    
    print(f"Total number of valid pairs: {valid_pairs_count}")

# Run the solver
find_coprime_prime_factor_pairs()
<<<0>>>
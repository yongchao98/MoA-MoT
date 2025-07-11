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
    Finds and verifies pairs (x,y) for the equation x^2-y^2=2023
    based on the given conditions.
    """
    N = 2023
    factor_pairs = []
    
    # Find all factor pairs (a,b) of N such that a > b > 0
    for b in range(1, int(math.sqrt(N)) + 1):
        if N % b == 0:
            a = N // b
            if a > b and b > 0:
                # Add pair (a,b) where a > b
                factor_pairs.append((a, b))
    
    # Reverse the list to check from largest 'a' to smallest
    factor_pairs.sort(reverse=True)
    
    print(f"The equation is x^2 - y^2 = 2023, which factors to (x+y)(x-y) = {N}.")
    print("Let a = x+y and b = x-y.")
    print(f"We need to find factor pairs (a, b) of {N} where a > b > 0, and check if they meet the conditions.\n")

    valid_pairs_count = 0

    for a, b in factor_pairs:
        print("----------------------------------------")
        print(f"Testing factor pair (a, b) = ({a}, {b})")
        
        # Step 1: Calculate x and y
        x = (a + b) // 2
        y = (a - b) // 2
        
        print("  Calculation of x and y:")
        print(f"    x = (a + b) / 2 = ({a} + {b}) / 2 = {x}")
        print(f"    y = (a - b) / 2 = ({a} - {b}) / 2 = {y}")

        # Step 2: Verify all conditions
        print("  Verification of conditions:")
        
        # Condition: x+y and x-y must be prime
        cond_a_is_prime = is_prime(a)
        cond_b_is_prime = is_prime(b)
        print(f"    1. Is x+y (= a = {a}) a prime number? {'Yes' if cond_a_is_prime else 'No'}")
        print(f"    2. Is x-y (= b = {b}) a prime number? {'Yes' if cond_b_is_prime else 'No'}")
        
        # Condition: x and y must be coprime
        gcd_xy = math.gcd(x, y)
        cond_coprime = (gcd_xy == 1)
        print(f"    3. Are x ({x}) and y ({y}) coprime? (gcd is {gcd_xy}) {'Yes' if cond_coprime else 'No'}")
        
        # Check if all conditions are met
        if cond_a_is_prime and cond_b_is_prime and cond_coprime:
            print(f"\n  Result: The pair (x, y) = ({x}, {y}) is a VALID solution.")
            valid_pairs_count += 1
        else:
            print(f"\n  Result: The pair (x, y) = ({x}, {y}) is INVALID as it fails one or more conditions.")

    print("----------------------------------------\n")
    print("Conclusion:")
    print(f"After checking all possible factor pairs, there are {valid_pairs_count} valid pairs (x,y) that satisfy all conditions.")
    print(f"<<<{valid_pairs_count}>>>")

# Run the main function
find_and_verify_pairs()
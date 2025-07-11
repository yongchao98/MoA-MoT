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
    Solves x^2 - y^2 = 2023 for coprime integers (x,y) where x+y and x-y are prime.
    """
    N = 2023
    print(f"Finding integer solutions for x^2 - y^2 = {N} where x > 0, y > 0.")
    print("This is equivalent to (x+y)(x-y) = N.")
    print("Let a = x+y and b = x-y. We need to find factors (a, b) of N such that a > b > 0.")
    print("Additional conditions: x and y are coprime, and both a and b must be prime numbers.\n")

    valid_pairs = []
    factor_pairs = []
    
    # Step 1: Find all factor pairs (a,b) of N with a > b
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            b_factor = i
            a_factor = N // i
            if a_factor > b_factor:
                factor_pairs.append((a_factor, b_factor))
    
    # The loop naturally finds pairs with descending 'a', so no sort is needed.

    print(f"Found {len(factor_pairs)} factor pairs (a, b) for {N} where a > b.")

    for i, (a, b) in enumerate(factor_pairs):
        print(f"\n--- Checking Factor Pair {i+1}: a = {a}, b = {b} ---")
        
        # Step 2: Verify if both a and b are prime
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        
        print(f"Condition 1: Check if a and b are prime.")
        print(f"Is a = {a} prime? {a_is_prime}")
        print(f"Is b = {b} prime? {b_is_prime}")

        if not (a_is_prime and b_is_prime):
            print("Result: This pair is not valid because both a and b must be prime.")
            continue

        # Step 3: Compute x and y
        # Since N is odd, its factors a and b must be odd.
        # So, a+b and a-b are even, which means x and y will be integers.
        x = (a + b) // 2
        y = (a - b) // 2

        print(f"\nCalculating x and y from the valid prime pair:")
        print(f"x = (a + b) / 2 = ({a} + {b}) / 2 = {x}")
        print(f"y = (a - b) / 2 = ({a} - {b}) / 2 = {y}")

        # Step 4: Verify if x and y are coprime
        are_coprime = (math.gcd(x, y) == 1)
        print(f"\nCondition 2: Check if x and y are coprime.")
        print(f"gcd({x}, {y}) = {math.gcd(x, y)}. Are they coprime? {are_coprime}")
        
        if are_coprime:
            print("\nResult: All conditions are met for this pair.")
            valid_pairs.append((x, y))
            print(f"Found a valid pair (x, y) = ({x}, {y}).")
            # Output the numbers in the final equation
            print(f"Final Equation: {x}^2 - {y}^2 = {N}")
        else:
            print("\nResult: This pair is not valid because x and y are not coprime.")

    print("\n------ FINAL SUMMARY ------")
    if not valid_pairs:
        print("No pairs (x,y) satisfy all the given conditions.")
    else:
        print("The following pairs (x,y) satisfy all conditions:")
        for pair in valid_pairs:
            print(pair)
            
    print(f"\nTotal number of valid pairs found: {len(valid_pairs)}")
    return len(valid_pairs)

if __name__ == "__main__":
    num_valid_pairs = find_and_verify_pairs()
    print(f"\n<<<{num_valid_pairs}>>>")
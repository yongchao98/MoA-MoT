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
    return a

def solve_equation():
    """
    Finds all pairs (x, y) satisfying the given conditions for x^2 - y^2 = 2023.
    """
    N = 2023
    print(f"Solving the equation x^2 - y^2 = {N}")
    print("The equation can be factored as (x+y)(x-y) = 2023.")
    print("Let a = x+y and b = x-y. We need to find factors (a, b) of 2023 where a > b > 0.\n")
    print("The conditions to check are:")
    print("1. Both a and b must be prime numbers.")
    print("2. x and y must be coprime integers.\n")

    # Step 1: Find all factor pairs of N where a > b
    factor_pairs = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            j = N // i
            # Ensure a > b
            if j > i:
                factor_pairs.append((j, i))
    
    print(f"Found {len(factor_pairs)} factor pair(s) (a, b) for {N} where a > b:")
    for a, b in factor_pairs:
        print(f"  - ({a}, {b})")
    print("\n--- Verifying conditions for each pair ---")

    valid_pairs_count = 0

    # Step 2: Iterate through factor pairs and check all conditions
    for i, (a, b) in enumerate(factor_pairs):
        print(f"\nCase {i+1}: Testing factor pair a = {a}, b = {b}")

        # Condition 1: a and b must be prime
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        print(f"Checking primality:")
        print(f"  - Is a = {a} prime? {a_is_prime}")
        print(f"  - Is b = {b} prime? {b_is_prime}")

        if not (a_is_prime and b_is_prime):
            print("Result: This pair is invalid because both 'a' and 'b' must be prime.\n")
            continue

        # This part of the code will only be reached if both a and b are prime.
        # Step 3: Calculate x and y
        x = (a + b) // 2
        y = (a - b) // 2
        print(f"Calculating x and y:")
        print(f"  - x = (a+b)/2 = ({a}+{b})/2 = {x}")
        print(f"  - y = (a-b)/2 = ({a}-{b})/2 = {y}")

        # Step 4: Check if x and y are coprime
        are_coprime = (gcd(x, y) == 1)
        print(f"Checking if x and y are coprime:")
        print(f"  - gcd({x}, {y}) = {gcd(x,y)}. Are they coprime? {are_coprime}")

        if are_coprime:
            print(f"Result: All conditions are met for (x, y) = ({x}, {y}). This is a valid pair.")
            valid_pairs_count += 1
            # Print the final equation as requested
            print(f"Verification of the equation: {x}^2 - {y}^2 = {x**2} - {y**2} = {x**2 - y**2}")
        else:
            print("Result: x and y are not coprime. This pair is invalid.")

    # Final summary
    print("\n--- Final Conclusion ---")
    if valid_pairs_count == 0:
        print("No pairs (x, y) were found that satisfy all the given conditions.")
    else:
        print(f"Found {valid_pairs_count} valid pair(s) that meet all conditions.")

    print(f"\nTotal number of valid pairs: {valid_pairs_count}")

if __name__ == '__main__':
    solve_equation()
<<<0>>>
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

def find_coprime_prime_factor_pairs():
    """
    Finds pairs (x, y) for x^2 - y^2 = N with special conditions.
    """
    N = 2023
    valid_pairs_count = 0

    print(f"The equation is x^2 - y^2 = {N}, which factors to (x + y)(x - y) = {N}.")
    print("Let a = x + y and b = x - y. We will find factor pairs (a, b) of 2023.")
    print("The conditions are: x, y are coprime, and both a and b must be prime numbers.\n")
    
    # Find all factor pairs (a, b) of N where a > b
    factor_pairs = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            # We need a > b, and since a*b=N, a will be N/i and b will be i.
            if (N // i) > i:
                factor_pairs.append((N // i, i))

    print(f"Found {len(factor_pairs)} factor pair(s) (a, b) of {N} where a > b:")
    for a, b in factor_pairs:
        print(f"  - ({a}, {b})")
    print("-" * 40)

    # Process each factor pair
    for a, b in factor_pairs:
        print(f"Processing pair (a, b) = ({a}, {b})")
        
        # Calculate x and y. a and b are both odd since their product is odd,
        # so x and y will be integers.
        x = (a + b) // 2
        y = (a - b) // 2
        
        print(f"  - Calculating x = (a+b)/2 = ({a}+{b})/2 = {x}")
        print(f"  - Calculating y = (a-b)/2 = ({a}-{b})/2 = {y}")

        # Verify all conditions
        is_a_prime = is_prime(a)
        is_b_prime = is_prime(b)
        are_coprime = math.gcd(x, y) == 1

        print("  - Verifying conditions:")
        print(f"    1. Is a = {a} a prime number? {'Yes.' if is_a_prime else 'No.'}")
        print(f"    2. Is b = {b} a prime number? {'Yes.' if is_b_prime else 'No.'}")
        print(f"    3. Are x = {x} and y = {y} coprime? {'Yes.' if are_coprime else 'No.'}")

        if is_a_prime and is_b_prime and are_coprime:
            valid_pairs_count += 1
            print("\n  - SUCCESS: All conditions met. This is a valid pair.")
            print(f"  - The final equation is: {x}^2 - {y}^2 = {x**2} - {y**2} = {N}")
        else:
            print("\n  - FAILURE: This pair is not valid because it fails one or more conditions.")
        
        print("-" * 40)

    print(f"Summary: After checking all factor pairs, the total number of valid pairs (x, y) is {valid_pairs_count}.")

# Run the main function
find_coprime_prime_factor_pairs()
<<<0>>>
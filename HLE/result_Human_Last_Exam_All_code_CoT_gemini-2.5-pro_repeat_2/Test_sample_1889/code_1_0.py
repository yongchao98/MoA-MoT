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
    while b:
        a, b = b, a % b
    return a

def find_and_verify_pairs():
    """
    Finds and verifies pairs (x,y) for the equation x^2 - y^2 = 2023
    based on the given conditions.
    """
    n = 2023
    valid_pairs_count = 0

    print(f"Solving x^2 - y^2 = {n} for coprime integers (x, y) where x, y > 0.")
    print("The equation can be factored as (x+y)(x-y) = 2023.")
    print("Let a = x+y and b = x-y. We are looking for factor pairs (a, b) of 2023 where a > b > 0.\n")

    # Find factor pairs (a, b) of n
    factor_pairs = []
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            # Ensure a > b
            a, b = n // i, i
            if a > b:
                factor_pairs.append((a, b))

    if not factor_pairs:
        print("No factor pairs found.")
    else:
        for a, b in factor_pairs:
            print(f"--- Checking factor pair (a={a}, b={b}) ---")
            
            # Calculate potential x and y
            x = (a + b) // 2
            y = (a - b) // 2
            print(f"  Calculating x = (a+b)/2 and y = (a-b)/2:")
            print(f"  x = ({a}+{b})/2 = {x}")
            print(f"  y = ({a}-{b})/2 = {y}\n")
            
            # Verify the conditions
            print(f"  Verifying conditions for the pair (x={x}, y={y}):")
            is_coprime_val = gcd(x, y) == 1
            a_is_prime = is_prime(a)
            b_is_prime = is_prime(b)

            print(f"    1. Are x and y coprime? gcd({x}, {y}) = {gcd(x, y)}.  -> {'Yes' if is_coprime_val else 'No'}")
            print(f"    2. Is x+y = {a} a prime number? -> {'Yes' if a_is_prime else 'No'}")
            print(f"    3. Is x-y = {b} a prime number? -> {'Yes' if b_is_prime else 'No'}")
            
            if is_coprime_val and a_is_prime and b_is_prime:
                valid_pairs_count += 1
                print("\n  RESULT: This is a valid pair.")
                print(f"  Final Equation: {x}^2 - {y}^2 = {x**2} - {y**2} = {n}")
            else:
                print("\n  RESULT: This pair is not valid as it fails one or more conditions.")
            print("-" * 45)

    print(f"\nSummary: After checking all factor pairs, the total number of valid pairs (x, y) is {valid_pairs_count}.")
    print(f"<<<{valid_pairs_count}>>>")

# Execute the function
find_and_verify_pairs()
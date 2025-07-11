import math

def is_prime(n):
    """
    Checks if a number n is prime.
    A prime number is a natural number greater than 1.
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

def find_solutions():
    """
    Finds and verifies solutions for the equation x^2 - y^2 = 2023
    based on the given conditions.
    """
    N = 2023
    valid_pairs_count = 0
    factor_pairs = []

    # Step 1: Find all pairs of factors (a, b) of N where a > b.
    # We let a = x+y and b = x-y.
    for b in range(1, int(math.sqrt(N)) + 1):
        if N % b == 0:
            a = N // b
            if a > b:
                factor_pairs.append((a, b))

    print(f"The equation is x^2 - y^2 = {N}")
    print(f"This is factored into (x+y)(x-y) = {N}.")
    print("We will examine all factor pairs (a, b) of 2023 where a = x+y and b = x-y.")
    
    # Iterate through factor pairs to find and verify solutions
    for a, b in sorted(factor_pairs, key=lambda p: p[1]):
        print(f"\n--- Checking Factor Pair (a={a}, b={b}) ---")
        
        # Step 2: Compute x and y
        x = (a + b) // 2
        y = (a - b) // 2
        print(f"Calculating potential solution (x,y):")
        print(f"  x = (a+b)/2 = ({a}+{b})/2 = {x}")
        print(f"  y = (a-b)/2 = ({a}-{b})/2 = {y}")
        
        # Step 3: Verify all conditions
        is_a_prime = is_prime(a)
        is_b_prime = is_prime(b)
        are_coprime = (math.gcd(x, y) == 1)

        print(f"Verifying conditions for (x={x}, y={y}):")
        print(f"  1. x+y and x-y must be prime:")
        print(f"     - Is x+y = {a} prime? {'Yes' if is_a_prime else 'No'}")
        print(f"     - Is x-y = {b} prime? {'Yes' if is_b_prime else 'No'}")
        print(f"  2. x and y must be coprime:")
        print(f"     - Is gcd({x}, {y}) = 1? {'Yes' if are_coprime else 'No'}")

        if is_a_prime and is_b_prime and are_coprime:
            valid_pairs_count += 1
            print("\nResult: SUCCESS! This pair satisfies all conditions.")
            # If a solution is found, print the equation with the numbers
            print(f"The valid equation is: {x}^2 - {y}^2 = {N}")
        else:
            print("\nResult: FAILURE. This pair does not satisfy all conditions.")

    print("\n--- Summary ---")
    print("How many valid pairs exist?")
    print(f"Based on the analysis, there are {valid_pairs_count} valid pairs (x,y).")

# Run the function
find_solutions()
print("<<<0>>>")
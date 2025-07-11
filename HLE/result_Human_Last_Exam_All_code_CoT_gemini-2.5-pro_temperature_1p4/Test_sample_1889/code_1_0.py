import math

def is_prime(n):
    """
    Checks if a number is prime.
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

def solve_equation():
    """
    Finds and verifies all integer pairs (x, y) for the equation x^2 - y^2 = 2023
    based on the specified conditions.
    """
    N = 2023
    print(f"Solving x^2 - y^2 = {N} for coprime integers x, y > 0.")
    print(f"This is equivalent to (x + y)(x - y) = {N}.")
    print("We let a = x+y and b = x-y, where a*b = N and a > b > 0.\n")

    valid_pairs_count = 0
    factor_pairs = []

    # Find all factor pairs (a, b) of N such that a > b
    for i in range(1, math.isqrt(N) + 1):
        if N % i == 0:
            b = i
            a = N // b
            if a > b:
                factor_pairs.append((a, b))

    # The problem asks to check from factors, so we process them systematically.
    # Sorting ensures a consistent order, from largest 'a' to smallest.
    factor_pairs.sort(key=lambda p: p[0], reverse=True)

    if not factor_pairs:
        print("No factor pairs (a,b) found for N={N} where a > b.")
    
    for i, (a, b) in enumerate(factor_pairs):
        print(f"--- Case {i+1}: Factoring 2023 as a = {a}, b = {b} ---")

        # Calculate x and y from a and b
        x = (a + b) // 2
        y = (a - b) // 2
        print(f"Calculating x and y:")
        print(f"x = (a + b) / 2 = ({a} + {b}) / 2 = {x}")
        print(f"y = (a - b) / 2 = ({a} - {b}) / 2 = {y}\n")
        
        print("Verifying the required conditions:")
        
        # 1. Check if x + y = a is prime
        is_a_prime = is_prime(a)
        print(f"1. Is x + y = {a} a prime number? {'Yes' if is_a_prime else 'No'}.")

        # 2. Check if x - y = b is prime
        is_b_prime = is_prime(b)
        print(f"2. Is x - y = {b} a prime number? {'Yes' if is_b_prime else 'No'}.")

        # 3. Check if x and y are coprime
        are_coprime = math.gcd(x, y) == 1
        print(f"3. Are x = {x} and y = {y} coprime? (gcd is {math.gcd(x, y)}) {'Yes' if are_coprime else 'No'}.")

        # Check if all conditions are met
        if is_a_prime and is_b_prime and are_coprime:
            valid_pairs_count += 1
            print(f"\nSUCCESS: All conditions met for (x, y) = ({x}, {y}).")
            # Output the equation as requested
            print(f"Equation check: {x}^2 - {y}^2 = {x**2} - {y**2} = {N}")
        else:
            print(f"\nFAILURE: The pair (x,y)=({x},{y}) is not a valid solution.")
        print("-" * 40 + "\n")

    print("-------------------- CONCLUSION --------------------")
    print(f"Total number of valid pairs found: {valid_pairs_count}")

# Run the solver
solve_equation()

# Final answer format
print("<<<0>>>")
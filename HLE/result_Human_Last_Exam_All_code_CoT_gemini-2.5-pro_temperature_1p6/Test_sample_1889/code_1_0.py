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

def find_solution_pairs():
    """
    Finds pairs (x,y) for x^2 - y^2 = N satisfying the given conditions.
    """
    N = 2023
    valid_pairs_count = 0
    
    print(f"Solving the equation x^2 - y^2 = {N} for coprime integers (x, y > 0).")
    print("This is equivalent to (x+y)(x-y) = 2023.")
    print("Let a = x+y and b = x-y. We need a and b to be prime numbers.\n")
    
    # Find factor pairs (a, b) of N where a > b
    factors = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            # Add the pair (a, b) where a > b
            factors.append((N // i, i))
            
    print(f"Prime factorization of {N} is 7 * 17^2.")
    print("The factor pairs (a, b) of 2023 are:")
    for a,b in factors:
        print(f"  - ({a}, {b})")
    print("-" * 30)

    # Check each factor pair against the conditions
    for a, b in factors:
        print(f"Checking factor pair a={a}, b={b}:")
        
        # Condition 1 & 2: a and b must be prime
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        
        print(f"  - Is 'a' ({a}) a prime number? {'Yes' if a_is_prime else 'No'}.")
        print(f"  - Is 'b' ({b}) a prime number? {'Yes' if b_is_prime else 'No'}.")
        
        if not (a_is_prime and b_is_prime):
            print("  - Result: Invalid pair, as both a and b must be prime.\n")
            continue
            
        # If we reached here, both a and b are prime.
        # Now we can calculate x and y and check for coprimality.
        x = (a + b) // 2
        y = (a - b) // 2
        
        # Condition 3: x and y must be coprime
        are_coprime = math.gcd(x, y) == 1
        
        print(f"  - Calculated (x,y) = ({x}, {y})")
        print(f"  - Are x and y coprime? {'Yes' if are_coprime else 'No'}.")
        
        if are_coprime:
            print("  - Result: This is a valid pair.")
            print(f"  - Final Equation: {x}^2 - {y}^2 = {N}")
            valid_pairs_count += 1
        else:
            print("  - Result: Invalid pair, as x and y are not coprime.\n")

    print("-" * 30)
    print(f"How many valid pairs exist?")
    print(f"The analysis found {valid_pairs_count} distinct pairs (x,y) satisfying all conditions.")

    return valid_pairs_count

# Run the process and get the final answer count
count = find_solution_pairs()
print(f"<<<{count}>>>")
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

def find_coprime_pairs_for_diophantine():
    """
    Solves x^2 - y^2 = 2023 for coprime positive integers (x, y)
    with the additional condition that x+y and x-y are prime.
    """
    N = 2023
    valid_pairs_count = 0
    
    print(f"Finding integer solutions for x^2 - y^2 = {N}")
    print("This can be factored as (x+y)(x-y) = 2023.")
    print("Let a = x+y and b = x-y.")
    print("We need to find factors (a, b) of 2023 such that a > b > 0.\n")
    
    # Find factor pairs of N where a > b
    factors = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            # Since N is odd, its factors i and N//i must also be odd.
            # a = N//i and b = i will always result in integer x and y.
            a = N // i
            b = i
            if a > b:
                factors.append((a, b))

    print(f"Found {len(factors)} factor pair(s) (a, b) for {N} where a > b.")
    
    # Process each factor pair
    for i, (a, b) in enumerate(factors):
        print(f"\n----- Case {i+1}: Factoring {N} = {a} * {b} -----")
        
        # Calculate x and y
        x = (a + b) // 2
        y = (a - b) // 2
        
        print(f"From a = x+y = {a} and b = x-y = {b}, we get:")
        print(f"x = (a+b)/2 = ({a}+{b})/2 = {x}")
        print(f"y = (a-b)/2 = ({a}-{b})/2 = {y}")
        
        # Verify conditions
        print("\nVerifying conditions:")
        
        # Condition 1: x+y and x-y are prime
        is_a_prime = is_prime(a)
        is_b_prime = is_prime(b)
        print(f"1. Are x+y and x-y prime?")
        print(f"   - Is x+y = {a} prime? {is_a_prime}")
        print(f"   - Is x-y = {b} prime? {is_b_prime}")
        
        # Condition 2: x and y are coprime
        are_coprime = (math.gcd(x, y) == 1)
        print(f"2. Are x = {x} and y = {y} coprime?")
        print(f"   - gcd({x}, {y}) = {math.gcd(x,y)}. Coprime? {are_coprime}")
        
        if is_a_prime and is_b_prime and are_coprime:
            print("\nResult: All conditions are satisfied. This is a valid pair.")
            print("The equation for this pair is:")
            print(f"{x}^2 - {y}^2 = {x*x} - {y*y} = {N}")
            valid_pairs_count += 1
        else:
            print("\nResult: Not all conditions are satisfied. This is not a valid pair.")

    print("\n-----------------------------------------")
    print(f"Total number of valid pairs (x,y) found: {valid_pairs_count}")

# Execute the function
find_coprime_pairs_for_diophantine()
print("<<<0>>>")
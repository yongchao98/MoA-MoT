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

def solve():
    """
    Finds pairs (x, y) for x^2 - y^2 = 2023 with the given conditions.
    """
    N = 2023
    valid_pairs_count = 0
    
    print(f"Finding integer solutions for x^2 - y^2 = {N}, where x > 0, y > 0.")
    print("The equation can be factored as (x+y)(x-y) = 2023.")
    print("Let a = x+y and b = x-y. We need to find factors (a,b) of 2023 where a > b > 0.\n")

    factor_pairs = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            factor1 = i
            factor2 = N // i
            # Ensure a > b
            a = max(factor1, factor2)
            b = min(factor1, factor2)
            if a != b:
                factor_pairs.append((a, b))

    print(f"Factor pairs (a, b) of {N} with a > b are: {factor_pairs}\n")

    for a, b in factor_pairs:
        print(f"--- Checking factor pair (a, b) = ({a}, {b}) ---")
        
        # Calculate x and y
        x = (a + b) // 2
        y = (a - b) // 2
        
        print(f"This yields x = ({a} + {b})/2 = {x} and y = ({a} - {b})/2 = {y}.")
        print(f"The corresponding equation is: {x}^2 - {y}^2 = {x*x} - {y*y} = {N}")

        # Verify the conditions
        print("\nVerifying conditions:")
        cond1_a_prime = is_prime(a)
        print(f"1. Is x+y = {a} a prime number? {cond1_a_prime}")
        
        cond2_b_prime = is_prime(b)
        print(f"2. Is x-y = {b} a prime number? {cond2_b_prime}")

        cond3_coprime = math.gcd(x, y) == 1
        print(f"3. Are x={x} and y={y} coprime (gcd=1)? {cond3_coprime}")

        if cond1_a_prime and cond2_b_prime and cond3_coprime:
            print("\n>>> SUCCESS: All conditions met. This is a valid pair.")
            valid_pairs_count += 1
        else:
            print("\n>>> FAILURE: This pair does not satisfy all conditions.")
        print("-" * 50)

    print(f"\nFinal Result: Found {valid_pairs_count} valid pair(s) of (x,y).")

solve()
print("<<<0>>>")
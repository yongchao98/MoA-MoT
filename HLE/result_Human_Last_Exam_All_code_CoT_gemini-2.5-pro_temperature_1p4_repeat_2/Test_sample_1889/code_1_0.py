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

def find_solution_pairs():
    """
    Solves x^2 - y^2 = 2023 with the given conditions.
    """
    N = 2023
    print(f"We are solving the equation x^2 - y^2 = {N}, where x>0, y>0.")
    print(f"This can be factored as (x+y)(x-y) = {N}.\n")
    print("Let a = x+y and b = x-y. We will find factors of {N} for a and b.")
    print("The conditions for a valid pair (x,y) are:")
    print("1. x and y are coprime (gcd(x, y) = 1).")
    print("2. x+y (which is 'a') must be a prime number.")
    print("3. x-y (which is 'b') must be a prime number.\n")

    # Find factor pairs (a, b) of N such that a > b
    factor_pairs = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            # We need a > b, so a = N//i and b = i
            if (N // i) > i:
                factor_pairs.append((N // i, i))

    print(f"Finding factor pairs (a, b) of {N} where a*b = {N} and a > b > 0:")
    
    valid_pairs_count = 0
    
    for a, b in factor_pairs:
        print(f"--- Checking factor pair (a={a}, b={b}) ---")
        
        # Calculate x and y from the factors a and b.
        # Since N=2023 is odd, its factors a and b must also be odd.
        # Thus (a+b) and (a-b) are guaranteed to be even, making x and y integers.
        x = (a + b) // 2
        y = (a - b) // 2
        
        print(f"From these factors, we calculate x and y:")
        print(f"x = ({a} + {b}) / 2 = {x}")
        print(f"y = ({a} - {b}) / 2 = {y}")

        # Check if all conditions are met
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        are_coprime = gcd(x, y) == 1
        
        print("\nVerifying the conditions:")
        print(f"1. Is x+y ({a}) a prime number? {'Yes' if a_is_prime else 'No'}")
        print(f"2. Is x-y ({b}) a prime number? {'Yes' if b_is_prime else 'No'}")
        print(f"3. Are x ({x}) and y ({y}) coprime? {'Yes' if are_coprime else 'No'} (gcd is {gcd(x,y)})")

        if a_is_prime and b_is_prime and are_coprime:
            valid_pairs_count += 1
            print(f"\nResult: The pair (x={x}, y={y}) is a VALID solution.")
        else:
            print(f"\nResult: The pair (x={x}, y={y}) is NOT a valid solution as it fails one or more conditions.")
        print("-" * 45 + "\n")

    print(f"Final Conclusion:")
    print(f"After checking all factor pairs, the total number of valid pairs (x,y) that satisfy all conditions is {valid_pairs_count}.")

if __name__ == '__main__':
    find_solution_pairs()
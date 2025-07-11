import math

def is_prime(n):
    """Checks if a number n is prime."""
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

def find_coprime_integer_pairs():
    """
    Finds all distinct pairs of coprime integers (x,y) for x^2-y^2=2023
    where x+y and x-y are prime.
    """
    N = 2023
    valid_pairs_count = 0

    print(f"Solving the equation x^2 - y^2 = {N}")
    print("This can be factored as (x+y)(x-y) = 2023.")
    print("Let a = x+y and b = x-y. We will check all factor pairs of 2023.")
    print("-" * 60)

    # Find all factor pairs (a, b) of N such that a > b > 0
    factor_pairs = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            a = N // i
            b = i
            if a > b:
                factor_pairs.append((a, b))

    # Check each factor pair against the conditions
    for a, b in factor_pairs:
        print(f"Checking factor pair (a, b) = ({a}, {b})")

        # Condition 1: Both a and b must be prime numbers.
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)

        if not (a_is_prime and b_is_prime):
            print(f"  - FAILED: Both numbers must be prime. (a is prime: {a_is_prime}, b is prime: {b_is_prime})")
            print("-" * 60)
            continue
        
        print(f"  - PASSED: Both a={a} and b={b} are prime.")

        # Calculate x and y. They are guaranteed to be integers since N is odd,
        # making its factors a and b odd, so their sum/difference is even.
        x = (a + b) // 2
        y = (a - b) // 2
        
        print(f"  - Calculating x = (a+b)/2 = {x}, y = (a-b)/2 = {y}")

        # Condition 2: x and y must be coprime.
        if gcd(x, y) == 1:
            print(f"  - PASSED: x={x} and y={y} are coprime (gcd=1).")
            print(f"  - SUCCESS! A valid pair (x,y) = ({x},{y}) has been found.")
            # Print the final equation with the numbers
            print(f"  - Verification: {x}^2 - {y}^2 = {x**2} - {y**2} = {x**2 - y**2}")
            valid_pairs_count += 1
        else:
            print(f"  - FAILED: x={x} and y={y} are not coprime (gcd={gcd(x, y)}).")
        
        print("-" * 60)

    print(f"\nConclusion: The total number of valid pairs is {valid_pairs_count}.")

if __name__ == '__main__':
    find_coprime_integer_pairs()
    print("<<<0>>>")
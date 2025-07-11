import math

def is_prime(n):
    """
    Checks if a number n is prime.
    A number is prime if it is greater than 1 and has no positive divisors other than 1 and itself.
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

def find_integer_pairs():
    """
    Finds all distinct pairs of coprime integers (x, y) for the equation x^2 - y^2 = 2023,
    with the additional constraints that x+y and x-y are prime numbers.
    """
    N = 2023
    print(f"We are solving the equation x^2 - y^2 = {N} for positive integers x and y.")
    print("This can be factored as (x + y)(x - y) = 2023.")
    print("Let a = x + y and b = x - y.")
    print("\nThe problem requires us to find pairs (x, y) that satisfy the following conditions:")
    print("1. x and y are positive integers.")
    print("2. x and y are coprime (their greatest common divisor is 1).")
    print("3. a = x + y is a prime number.")
    print("4. b = x - y is a prime number.")
    print(f"\nFirst, we need to find pairs of factors (a, b) of {N} such that a * b = {N} and a > b.")

    # Find all pairs of factors of N
    factor_pairs = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            # We have a pair of factors (i, N/i)
            # Since x,y > 0, we must have x+y > x-y > 0. So a > b > 0.
            b = i
            a = N // i
            if a > b:
                factor_pairs.append((a, b))

    print(f"\nThe pairs of factors (a, b) of {N} where a > b are: {factor_pairs}\n")

    valid_pairs_count = 0

    # For each factor pair, check if they meet the conditions
    for a, b in factor_pairs:
        print(f"--- Verifying Factor Pair (a = {a}, b = {b}) ---")
        
        # Condition check: a and b must both be prime
        print(f"Step 1: Check if both a={a} and b={b} are prime numbers.")
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)

        if not a_is_prime:
            print(f"Result: a = {a} is NOT prime. (2023 = 7 * 17^2)")
        else:
            print(f"Result: a = {a} IS prime.")
        
        if not b_is_prime:
            print(f"Result: b = {b} is NOT prime.")
        else:
            print(f"Result: b = {b} IS prime.")

        if not (a_is_prime and b_is_prime):
            print("Conclusion: Since at least one of a or b is not prime, this pair is invalid.\n")
            continue
        
        # The following block of code would execute if a pair of prime factors were found.
        # Based on the factors of 2023, this part will not be reached.
        print("\nStep 2: Since both a and b are prime, we can calculate x and y.")
        x = (a + b) // 2
        y = (a - b) // 2
        
        print(f"x = (a + b) / 2 = ({a} + {b}) / 2 = {x}")
        print(f"y = (a - b) / 2 = ({a} - {b}) / 2 = {y}")

        # Output the numbers in the final equation as requested
        print(f"Proposed equation: {x}^2 - {y}^2 = ({x}+{y})*({x}-{y}) = {a}*{b} = {N}")

        print("\nStep 3: Check if x and y are coprime.")
        if math.gcd(x, y) == 1:
            print(f"Result: gcd({x}, {y}) is 1. They are coprime.")
            print("Conclusion: All conditions are satisfied. This is a valid pair.")
            valid_pairs_count += 1
        else:
            print(f"Result: gcd({x}, {y}) is not 1. They are NOT coprime.")
            print("Conclusion: This pair is invalid.\n")
    
    print("--- Final Result ---")
    print(f"After checking all factor pairs of {N}, we found {valid_pairs_count} valid pairs (x, y) that satisfy all the given conditions.")
    print(f"\nHow many valid pairs exist? Answer: {valid_pairs_count}")

if __name__ == '__main__':
    find_integer_pairs()
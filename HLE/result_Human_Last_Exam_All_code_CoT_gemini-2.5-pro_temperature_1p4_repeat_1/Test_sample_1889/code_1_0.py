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

def find_coprime_pairs():
    """
    Finds all distinct pairs of coprime integers (x, y) for the equation x^2 - y^2 = 2023,
    with the additional condition that x+y and x-y must be prime numbers.
    """
    N = 2023
    valid_pairs_count = 0

    print(f"Solving the equation x^2 - y^2 = {N}")
    print("This can be factored into (x+y)(x-y) = 2023.")
    print("Let a = x+y and b = x-y. We are looking for factors of 2023.")
    print("The given conditions require that x and y are coprime, and both a and b are prime numbers.\n")

    # Step 1: Find all pairs of factors (a, b) for N such that a > b > 0.
    factor_pairs = []
    # We only need to check up to the square root of N.
    for i in range(1, math.isqrt(N) + 1):
        if N % i == 0:
            b = i
            a = N // i
            # We require a > b, so we only add this pair.
            if a > b:
                factor_pairs.append((a, b))

    print(f"The pairs of integer factors (a, b) for {N} where a > b are: {factor_pairs}\n")

    # Step 2: For each factor pair, check if they meet the conditions.
    for a, b in factor_pairs:
        print(f"--- Checking factor pair (a={a}, b={b}) ---")
        
        # Condition 1: a and b must be prime numbers.
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        
        print(f"Is a = {a} a prime number? {a_is_prime}")
        print(f"Is b = {b} a prime number? {b_is_prime}")

        if not (a_is_prime and b_is_prime):
            print("Conclusion: At least one of a or b is not prime. This pair is invalid.\n")
            continue
            
        # The following steps are for pairs where both a and b are prime.
        # Based on the factorization of 2023, this section will not be reached.
        
        # Step 3: Calculate x and y.
        x = (a + b) // 2
        y = (a - b) // 2
        
        # Step 4: Check if x and y are coprime.
        # If a and b are distinct primes, gcd(a,b)=1, which implies gcd(x,y)=1.
        if math.gcd(x, y) == 1:
            print(f"Valid pair found: (x, y) = ({x}, {y})")
            print(f"Verification: x+y = {x+y} (prime), x-y = {x-y} (prime), gcd(x,y) = 1.")
            # Print the final equation with the numbers
            print(f"Equation: {x}^2 - {y}^2 = {x*x} - {y*y} = {N}")
            valid_pairs_count += 1
        else:
            # This case won't be reached if a and b are distinct primes.
            print(f"x={x} and y={y} are not coprime. This pair is invalid.\n")

    print("--- Final Result ---")
    print(f"The total number of valid pairs (x, y) that satisfy all conditions is {valid_pairs_count}.")

if __name__ == '__main__':
    find_coprime_pairs()
    
<<<0>>>
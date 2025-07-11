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

def find_factor_pairs(n):
    """Finds all pairs of factors (a, b) of n such that a > b."""
    pairs = []
    # Iterate up to the square root of n to find factors
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            b = i
            a = n // i
            # Ensure a > b to avoid duplicate pairs and satisfy x > y
            if a > b:
                pairs.append((a, b))
    # Reverse to have pairs with larger 'a' first
    return sorted(pairs, key=lambda p: p[0], reverse=True)

def main():
    """
    Finds pairs of coprime integers (x, y) for x^2 - y^2 = 2023
    with the additional condition that x+y and x-y are prime.
    """
    N = 2023
    print(f"The equation is x^2 - y^2 = {N}")
    print(f"Factoring gives (x + y)(x - y) = {N}.")
    print("Let a = x + y and b = x - y. We will find factor pairs (a,b) of 2023.\n")

    factor_pairs = find_factor_pairs(N)
    valid_pairs_count = 0

    print(f"The factor pairs (a,b) of {N} with a > b are: {factor_pairs}\n")

    for a, b in factor_pairs:
        print(f"--- Checking factor pair (a={a}, b={b}) ---")

        # Calculate x and y. a,b are odd factors of an odd number,
        # so a+b and a-b will be even, and x, y will be integers.
        x = (a + b) // 2
        y = (a - b) // 2

        print(f"Calculated (x,y): x = ({a}+{b})/2 = {x}, y = ({a}-{b})/2 = {y}")
        
        # Outputting the numbers in the equation
        print(f"Checking the equation with this pair: {x}^2 - {y}^2 = {x**2 - y**2}")
        print(f"Factored form: ({x}+{y}) * ({x}-{y}) = {x+y} * {x-y} = {N}")


        print("\nVerifying all conditions:")
        cond_a_prime = is_prime(a)
        cond_b_prime = is_prime(b)
        # It can be proven that GCD(x,y) = GCD(a,b) when a,b are odd.
        coprime = math.gcd(x, y) == 1

        print(f"1. x+y = {a} is prime? {'Yes' if cond_a_prime else 'No'}")
        print(f"2. x-y = {b} is prime? {'Yes' if cond_b_prime else 'No'}")
        print(f"3. x ({x}) and y ({y}) are coprime? {'Yes' if coprime else 'No'} (GCD={math.gcd(x,y)})")

        if cond_a_prime and cond_b_prime and coprime:
            valid_pairs_count += 1
            print("\nResult: VALID PAIR FOUND")
        else:
            print("\nResult: This pair is invalid.")
        print("-" * 40 + "\n")

    print(f"Conclusion: After checking all possibilities, we found {valid_pairs_count} valid pairs.")

if __name__ == '__main__':
    main()
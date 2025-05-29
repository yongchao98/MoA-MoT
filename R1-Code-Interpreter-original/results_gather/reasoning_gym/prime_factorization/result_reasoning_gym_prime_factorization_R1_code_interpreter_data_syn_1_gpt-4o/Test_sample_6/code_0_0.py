def prime_factors(n):
    factors = []
    # Start with the smallest prime number
    divisor = 2
    while n >= 2:
        if n % divisor == 0:
            factors.append(divisor)
            n = n // divisor
        else:
            divisor += 1
    return factors

# Find the prime factors of 656
factors = prime_factors(656)
print(" × ".join(map(str, factors)))
def prime_factors(n):
    factors = []
    # Start with the smallest prime number
    divisor = 2
    while n > 1:
        while n % divisor == 0:
            factors.append(divisor)
            n //= divisor
        divisor += 1
    return factors

factors = prime_factors(656)
print(" × ".join(map(str, factors)))
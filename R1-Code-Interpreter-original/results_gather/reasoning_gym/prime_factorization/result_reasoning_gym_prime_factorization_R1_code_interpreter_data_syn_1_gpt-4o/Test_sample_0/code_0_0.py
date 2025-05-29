def prime_factors(n):
    factors = []
    # Start with the smallest prime number
    divisor = 2
    while n >= divisor:
        while n % divisor == 0:
            factors.append(divisor)
            n //= divisor
        divisor += 1
    return factors

# Find the prime factors of 81
factors_of_81 = prime_factors(81)
# Format the output as required
formatted_factors = ' × '.join(map(str, factors_of_81))
print(formatted_factors)
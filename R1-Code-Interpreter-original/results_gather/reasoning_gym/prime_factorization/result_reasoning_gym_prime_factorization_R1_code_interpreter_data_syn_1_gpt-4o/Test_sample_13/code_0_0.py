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

# Find the prime factors of 280
factors_of_280 = prime_factors(280)
# Format the output as required
formatted_factors = ' × '.join(map(str, factors_of_280))
print(f"<<<{formatted_factors}>>>")
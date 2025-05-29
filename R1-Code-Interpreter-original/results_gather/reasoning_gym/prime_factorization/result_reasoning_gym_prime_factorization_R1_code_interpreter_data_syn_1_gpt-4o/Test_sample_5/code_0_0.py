def prime_factors(n):
    factors = []
    # Start with the smallest prime number
    divisor = 2
    while n >= 2:
        while n % divisor == 0:
            factors.append(divisor)
            n //= divisor
        divisor += 1
    return factors

# Find the prime factors of 700
factors_of_700 = prime_factors(700)
# Format the output as required
formatted_factors = ' × '.join(map(str, factors_of_700))
print(f"<<<{formatted_factors}>>>")
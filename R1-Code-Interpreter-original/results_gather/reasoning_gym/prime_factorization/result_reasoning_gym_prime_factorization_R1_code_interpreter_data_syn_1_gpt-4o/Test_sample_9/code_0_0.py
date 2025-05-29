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

# Find the prime factors of 392
factors = prime_factors(392)
# Format the output as requested
result = " × ".join(map(str, factors))
print(f"<<<{result}>>>")
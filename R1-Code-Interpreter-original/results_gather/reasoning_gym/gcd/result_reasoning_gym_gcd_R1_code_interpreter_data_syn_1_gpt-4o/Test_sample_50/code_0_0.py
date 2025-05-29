def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Find the GCD of 788 and 502
result = gcd(788, 502)
print(result)
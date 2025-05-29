def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Find the GCD of 24 and 140
result = gcd(24, 140)
print(result)
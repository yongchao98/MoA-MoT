def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Find the GCD of 56 and 902
result = gcd(56, 902)
print(result)
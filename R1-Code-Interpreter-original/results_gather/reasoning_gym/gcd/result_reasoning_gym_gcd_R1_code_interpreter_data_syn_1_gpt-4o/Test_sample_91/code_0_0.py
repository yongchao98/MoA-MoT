def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Find the GCD of 275 and 90
result = gcd(275, 90)
print(result)
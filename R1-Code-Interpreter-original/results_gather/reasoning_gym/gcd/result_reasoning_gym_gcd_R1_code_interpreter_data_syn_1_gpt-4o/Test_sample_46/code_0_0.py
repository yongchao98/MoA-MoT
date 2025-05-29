def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Find the GCD of 99 and 308
result = gcd(99, 308)
print(result)
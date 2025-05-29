def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Find the GCD of 642 and 640
result = gcd(642, 640)
print(result)
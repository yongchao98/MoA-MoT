def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Find the GCD of 302 and 720
result = gcd(302, 720)
print(result)
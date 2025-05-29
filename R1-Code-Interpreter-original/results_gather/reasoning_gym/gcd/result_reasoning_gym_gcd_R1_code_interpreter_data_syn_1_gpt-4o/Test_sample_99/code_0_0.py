def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Find the GCD of 400 and 798
result = gcd(400, 798)
print(result)
def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Find the GCD of 758 and 310
result = gcd(758, 310)
print(result)
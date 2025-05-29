def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Find the GCD of 159 and 68
result = gcd(159, 68)
print(result)
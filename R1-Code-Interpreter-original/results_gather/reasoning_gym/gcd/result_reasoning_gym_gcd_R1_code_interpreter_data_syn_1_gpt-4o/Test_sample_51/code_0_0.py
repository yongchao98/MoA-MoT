def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Find the GCD of 36 and 51
result = gcd(36, 51)
print(result)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Find the GCD of 60 and 22
gcd_value = gcd(60, 22)
print(gcd_value)
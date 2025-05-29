def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Find the GCD of 26 and 760
gcd_value = gcd(26, 760)
print(gcd_value)
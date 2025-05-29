def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Find the GCD of 944 and 966
result = gcd(944, 966)
print(result)
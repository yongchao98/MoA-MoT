def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Numbers to find the GCD of
num1 = 55
num2 = 429

# Calculate GCD
result = gcd(num1, num2)
print(result)
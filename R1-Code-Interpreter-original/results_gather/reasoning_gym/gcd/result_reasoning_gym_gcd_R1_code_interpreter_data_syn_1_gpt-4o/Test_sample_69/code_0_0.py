def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Given numbers
num1 = 617
num2 = 228

# Calculate GCD
result = gcd(num1, num2)
print(result)
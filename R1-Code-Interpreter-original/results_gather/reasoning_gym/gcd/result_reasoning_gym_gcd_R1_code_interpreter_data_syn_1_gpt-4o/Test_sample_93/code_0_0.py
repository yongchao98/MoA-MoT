def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Given numbers
num1 = 228
num2 = 849

# Calculate GCD
result = gcd(num1, num2)
print(result)
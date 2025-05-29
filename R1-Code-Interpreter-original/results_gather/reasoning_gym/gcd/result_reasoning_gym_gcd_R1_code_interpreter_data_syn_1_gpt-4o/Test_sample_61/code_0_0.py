def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Given numbers
num1 = 364
num2 = 152

# Calculate GCD
result = gcd(num1, num2)
print(result)
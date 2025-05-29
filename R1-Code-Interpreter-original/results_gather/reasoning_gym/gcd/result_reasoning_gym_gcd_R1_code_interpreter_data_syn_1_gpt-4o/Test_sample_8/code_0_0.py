def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Given numbers
num1 = 252
num2 = 290

# Calculate GCD
result = gcd(num1, num2)
print(result)
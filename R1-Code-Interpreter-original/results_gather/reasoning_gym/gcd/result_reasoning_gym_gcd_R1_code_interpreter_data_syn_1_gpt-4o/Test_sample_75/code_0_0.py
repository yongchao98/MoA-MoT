def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Given numbers
num1 = 980
num2 = 142

# Calculate GCD
result = gcd(num1, num2)
print(result)
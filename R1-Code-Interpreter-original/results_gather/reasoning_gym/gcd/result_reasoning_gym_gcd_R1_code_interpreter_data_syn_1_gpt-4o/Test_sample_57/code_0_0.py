def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Given numbers
num1 = 154
num2 = 868

# Calculate GCD
result = gcd(num1, num2)
print(result)
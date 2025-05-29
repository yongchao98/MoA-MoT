def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Given numbers
num1 = 949
num2 = 219

# Calculate GCD
result = gcd(num1, num2)
print(result)
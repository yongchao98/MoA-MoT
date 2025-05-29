def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Given numbers
num1 = 450
num2 = 570

# Calculate GCD
result = gcd(num1, num2)
print(result)
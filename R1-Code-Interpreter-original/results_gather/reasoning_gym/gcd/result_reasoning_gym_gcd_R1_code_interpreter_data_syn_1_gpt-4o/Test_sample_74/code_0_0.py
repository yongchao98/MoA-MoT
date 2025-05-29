def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Given numbers
num1 = 824
num2 = 66

# Calculate GCD
result = gcd(num1, num2)
print(result)
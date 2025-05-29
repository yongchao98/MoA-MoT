def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Given numbers
num1 = 274
num2 = 968

# Calculate GCD
result = gcd(num1, num2)
print(result)
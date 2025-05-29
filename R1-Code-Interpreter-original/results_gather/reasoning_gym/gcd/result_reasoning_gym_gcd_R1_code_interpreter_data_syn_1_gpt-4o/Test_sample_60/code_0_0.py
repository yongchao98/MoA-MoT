def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Given numbers
num1 = 584
num2 = 536

# Calculate GCD
result = gcd(num1, num2)
print(result)
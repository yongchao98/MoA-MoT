def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Given numbers
num1 = 732
num2 = 781

# Calculate GCD
result = gcd(num1, num2)
print(result)
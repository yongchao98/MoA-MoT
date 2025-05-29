def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Given numbers
num1 = 390
num2 = 205

# Calculate GCD
gcd_result = gcd(num1, num2)
print(gcd_result)
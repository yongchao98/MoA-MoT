def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Given numbers
num1 = 198
num2 = 843

# Calculate GCD
gcd_result = gcd(num1, num2)
print(gcd_result)
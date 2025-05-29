def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Given numbers
num1 = 697
num2 = 46

# Calculate GCD
gcd_result = gcd(num1, num2)
print(gcd_result)
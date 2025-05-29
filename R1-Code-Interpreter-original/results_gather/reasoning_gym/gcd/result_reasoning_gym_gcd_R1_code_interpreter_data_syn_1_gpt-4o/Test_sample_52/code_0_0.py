def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Given numbers
num1 = 936
num2 = 794

# Calculate GCD
result = gcd(num1, num2)
print(result)
def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Given numbers
num1 = 853
num2 = 566

# Calculate GCD
result = gcd(num1, num2)
print(result)
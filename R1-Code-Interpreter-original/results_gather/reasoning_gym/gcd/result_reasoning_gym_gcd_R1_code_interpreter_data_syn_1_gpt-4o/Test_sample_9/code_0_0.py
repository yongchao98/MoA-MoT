def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Numbers to find the GCD of
num1 = 924
num2 = 936

# Calculate GCD
result = gcd(num1, num2)
print(result)
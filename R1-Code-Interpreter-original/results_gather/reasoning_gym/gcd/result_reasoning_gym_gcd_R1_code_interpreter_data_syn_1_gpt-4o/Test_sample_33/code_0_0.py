def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# Numbers to find the GCD of
num1 = 455
num2 = 470

# Calculate GCD
gcd_result = gcd(num1, num2)
print(gcd_result)
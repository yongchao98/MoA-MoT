def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

result = gcd(320, 102)
print(result)
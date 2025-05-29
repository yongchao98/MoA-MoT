def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

result = gcd(200, 226)
print(result)
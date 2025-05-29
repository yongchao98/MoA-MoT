import math

numerator = 7476
denominator = 8388

gcd = math.gcd(numerator, denominator)

simplified_numerator = numerator // gcd
simplified_denominator = denominator // gcd

print(f"{simplified_numerator}/{simplified_denominator}")
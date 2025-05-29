import math

numerator = 1978
denominator = 13915

gcd = math.gcd(numerator, denominator)

simplified_numerator = numerator // gcd
simplified_denominator = denominator // gcd

print(f"{simplified_numerator}/{simplified_denominator}")
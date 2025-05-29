import math

numerator = 398
denominator = 957

gcd = math.gcd(numerator, denominator)

simplified_numerator = numerator // gcd
simplified_denominator = denominator // gcd

print(f"{simplified_numerator}/{simplified_denominator}")
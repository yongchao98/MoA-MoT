import math

numerator = 3188
denominator = 3244

gcd = math.gcd(numerator, denominator)

simplified_numerator = numerator // gcd
simplified_denominator = denominator // gcd

print(f"{simplified_numerator}/{simplified_denominator}")
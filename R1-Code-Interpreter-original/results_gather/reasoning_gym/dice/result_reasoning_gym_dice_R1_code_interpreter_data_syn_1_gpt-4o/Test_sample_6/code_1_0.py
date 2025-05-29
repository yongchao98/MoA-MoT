import math

numerator = 17635
denominator = 22500

gcd = math.gcd(numerator, denominator)
simplified_numerator = numerator // gcd
simplified_denominator = denominator // gcd

print(f"{simplified_numerator}/{simplified_denominator}")
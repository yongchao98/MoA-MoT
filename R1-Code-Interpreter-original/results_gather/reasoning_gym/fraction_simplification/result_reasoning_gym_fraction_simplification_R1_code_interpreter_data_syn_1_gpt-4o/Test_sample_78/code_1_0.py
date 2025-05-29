# Check if the fraction 398/957 can be simplified further
numerator = 398
denominator = 957

gcd = math.gcd(numerator, denominator)

simplified_numerator = numerator // gcd
simplified_denominator = denominator // gcd

print(f"{simplified_numerator}/{simplified_denominator}")
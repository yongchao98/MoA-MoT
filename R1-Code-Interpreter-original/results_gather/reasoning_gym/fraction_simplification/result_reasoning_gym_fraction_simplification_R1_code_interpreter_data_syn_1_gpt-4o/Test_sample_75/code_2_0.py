import math

# Initial simplified fraction from previous step
simplified_numerator = 253
simplified_denominator = 635

# Further simplify the fraction
gcd_final = math.gcd(simplified_numerator, simplified_denominator)

# Simplify the fraction further
final_numerator = simplified_numerator // gcd_final
final_denominator = simplified_denominator // gcd_final

# Output the fully simplified fraction
print(f"{final_numerator}/{final_denominator}")
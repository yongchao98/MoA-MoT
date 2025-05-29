import math

# Define the simplified numerator and denominator from the previous result
simplified_numerator = 65
simplified_denominator = 69

# Calculate the GCD of the simplified fraction
gcd_simplified = math.gcd(simplified_numerator, simplified_denominator)

# Further simplify the fraction if possible
final_numerator = simplified_numerator // gcd_simplified
final_denominator = simplified_denominator // gcd_simplified

# Print the final simplified fraction
print(f"{final_numerator}/{final_denominator}")
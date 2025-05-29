import math

# Define the numerator and denominator
numerator = 35970
denominator = 51414

# Calculate the GCD
gcd = math.gcd(numerator, denominator)

# Simplify the fraction
simplified_numerator = numerator // gcd
simplified_denominator = denominator // gcd

# Print the simplified fraction
print(f"{simplified_numerator}/{simplified_denominator}")
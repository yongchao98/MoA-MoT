import math

# The equation for the generating amplitude c1 (in the case c1=c2) is:
# 2*pi - (c1^2 / 2) * (e^(4*pi) - 1) = 0
# We solve for c1 > 0.

# Rearranging the equation to solve for c1^2:
# c1^2 = 4*pi / (e^(4*pi) - 1)

pi_val = math.pi
four_pi = 4 * pi_val

# Calculate the numerator and denominator of the fraction
numerator = four_pi
denominator = math.exp(four_pi) - 1

# Calculate c1^2
c1_squared = numerator / denominator

# Calculate c1 by taking the positive square root
c1 = math.sqrt(c1_squared)

print("The equation for the generating amplitude c1 (with c1=c2) is:")
print("2*pi - (c1^2 / 2) * (e^(4*pi) - 1) = 0")
print("\nSolving for c1, we get:")
print("c1 = sqrt(4*pi / (e^(4*pi) - 1))")
print("\nSubstituting the value of pi:")
print(f"c1 = sqrt((4 * {pi_val}) / (e^(4 * {pi_val}) - 1))")
print(f"c1 = sqrt({numerator} / ({math.exp(four_pi)} - 1))")
print(f"c1 = sqrt({numerator} / {denominator})")
print(f"c1 = sqrt({c1_squared})")
print(f"\nThe value of the first positive root c1 is:")
print(c1)
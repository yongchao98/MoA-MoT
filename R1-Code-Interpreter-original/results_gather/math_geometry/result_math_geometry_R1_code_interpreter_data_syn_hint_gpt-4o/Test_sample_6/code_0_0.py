import math

# Given values for Triangle 2
a = 3
b = 5
C = 120  # angle in degrees

# Convert angle to radians for cosine function
C_rad = math.radians(C)

# Law of Cosines to find c (AB)
AB = math.sqrt(a**2 + b**2 - 2 * a * b * math.cos(C_rad))

print(AB)
import math

# Given values
r = 3
h = 10

# Calculate lateral surface area
A_lateral = 2 * math.pi * r * h

# Calculate area of the two ends
A_ends = 2 * math.pi * r**2

# Calculate total surface area
A_total = A_lateral + A_ends

# Output the total surface area
print(A_total)
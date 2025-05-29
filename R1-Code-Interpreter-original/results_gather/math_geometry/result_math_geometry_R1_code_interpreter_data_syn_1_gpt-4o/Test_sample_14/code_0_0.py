import math

# Define the radii
radii = [2, 3, 5]

# Calculate the total volume
total_volume = sum((4/3) * math.pi * (r**3) for r in radii)

# Print the total volume in terms of pi
print(total_volume / math.pi)
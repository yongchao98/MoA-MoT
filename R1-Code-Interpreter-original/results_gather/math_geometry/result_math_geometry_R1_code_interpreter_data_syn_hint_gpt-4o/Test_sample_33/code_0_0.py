import math

# Given face areas
xy = 24
yz = 16
zx = 6

# Calculate the volume
volume_squared = xy * yz * zx
volume = math.sqrt(volume_squared)

print(volume)
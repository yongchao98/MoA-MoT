import math

# Given values
radius = 58
chord_length = 80

# Using the cosine rule to find the angle C
cos_C = (2 * radius**2 - chord_length**2) / (2 * radius**2)
angle_C = math.acos(cos_C)

# Using the angle to find the length of the second part of the journey (x)
# The second part of the journey is the arc length
x = 2 * radius * math.sin(angle_C / 2)

# Total distance traveled
total_distance = 116 + x + 80

print(total_distance)
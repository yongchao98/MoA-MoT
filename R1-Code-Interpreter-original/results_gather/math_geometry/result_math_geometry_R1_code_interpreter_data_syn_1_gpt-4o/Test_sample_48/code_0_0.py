import math

# Given values
diameter = 2
radius = diameter / 2
height_glass = 6
height_lemonade = height_glass / 2

# Volume of lemonade in the glass
volume_lemonade = math.pi * (radius ** 2) * height_lemonade

# Ratio of lemon juice to total parts
ratio_lemon_juice = 1 / 12

# Volume of lemon juice
volume_lemon_juice = volume_lemonade * ratio_lemon_juice

# Print the volume of lemon juice rounded to the nearest hundredth
print(round(volume_lemon_juice, 2))
import math

# Given values
radius = 1  # in inches
height_of_lemonade = 3  # in inches
ratio_lemon_juice = 1
ratio_water = 11

# Calculate the volume of lemonade
volume_of_lemonade = math.pi * (radius ** 2) * height_of_lemonade

# Calculate the volume of lemon juice
volume_of_lemon_juice = (ratio_lemon_juice / (ratio_lemon_juice + ratio_water)) * volume_of_lemonade

# Print the volume of lemon juice rounded to the nearest hundredth
print(round(volume_of_lemon_juice, 2))
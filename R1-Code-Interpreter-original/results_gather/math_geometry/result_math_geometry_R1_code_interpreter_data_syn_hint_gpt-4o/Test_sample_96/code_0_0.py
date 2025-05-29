import math

# Radii of the circles
r_large = 5.5
r_small1 = 2
r_small2 = 3

# Areas of the circles
area_large = math.pi * r_large**2
area_small1 = math.pi * r_small1**2
area_small2 = math.pi * r_small2**2

# Shaded area
shaded_area = area_large - area_small1 - area_small2

print(shaded_area)
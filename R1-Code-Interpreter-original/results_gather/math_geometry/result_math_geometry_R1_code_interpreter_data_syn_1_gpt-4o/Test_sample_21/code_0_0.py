# Constants
diameter = 2  # inches
radius = diameter / 2
length_in_inches = 12  # 1 foot = 12 inches

# Area of one semicircle
area_of_one_semicircle = (1/2) * 3.141592653589793 * (radius ** 2)

# Number of semicircles in 1 foot
number_of_semicircles = length_in_inches / diameter

# Total shaded area
total_shaded_area = number_of_semicircles * area_of_one_semicircle

print(total_shaded_area)
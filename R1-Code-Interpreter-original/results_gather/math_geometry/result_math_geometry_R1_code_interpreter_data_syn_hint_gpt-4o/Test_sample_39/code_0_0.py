# Given data
value_2_inch_cube = 200  # in dollars
volume_2_inch_cube = 2**3  # in cubic inches

# Calculate value per cubic inch
value_per_cubic_inch = value_2_inch_cube / volume_2_inch_cube

# Volume of the 3-inch cube
volume_3_inch_cube = 3**3  # in cubic inches

# Calculate the value of the 3-inch cube
value_3_inch_cube = value_per_cubic_inch * volume_3_inch_cube

# Round to the nearest dollar
value_3_inch_cube_rounded = round(value_3_inch_cube)

print(value_3_inch_cube_rounded)
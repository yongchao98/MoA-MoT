# Given data
side_2_inch = 2
weight_2_inch = 3  # pounds
value_2_inch = 200  # dollars

# Calculate volume of 2-inch cube
volume_2_inch = side_2_inch ** 3

# Calculate value per cubic inch
value_per_cubic_inch = value_2_inch / volume_2_inch

# Calculate volume of 3-inch cube
side_3_inch = 3
volume_3_inch = side_3_inch ** 3

# Calculate value of 3-inch cube
value_3_inch = value_per_cubic_inch * volume_3_inch

# Round to the nearest dollar
value_3_inch_rounded = round(value_3_inch)

print(value_3_inch_rounded)
# Python code to calculate the real cube root and round to 4 decimal places
import math

# Calculate the cube root of a negative number
x_cubed = -19 / 27
x = -math.pow(abs(x_cubed), 1/3)
rounded_x = round(x, 4)
print(rounded_x)
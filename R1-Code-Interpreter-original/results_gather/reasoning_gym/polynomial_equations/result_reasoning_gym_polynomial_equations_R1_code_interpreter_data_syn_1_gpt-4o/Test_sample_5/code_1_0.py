# Python code to calculate the real cube root and round to 4 decimal places
import cmath

# Calculate the cube root using cmath to handle negative numbers correctly
x_cubed = -19 / 27
x = cmath.exp(cmath.log(x_cubed) / 3).real
rounded_x = round(x, 4)
print(rounded_x)
import math

# Step 1: Define the parameters of the self-affine fractal.
# N is the number of recursive copies (the number of black keys).
N = 5

# r_y is the scaling factor in the height dimension (the larger scaling factor).
# A black key's height is 9/14 of the original keyboard's height.
r_y_num = 9
r_y_den = 14

# r_x is the scaling factor in the width dimension (the smaller scaling factor).
# A black key's width is (3/14), and the original keyboard's width is 3.
# So the scaling is (3/14) / 3 = 1/14.
r_x_num = 1
r_x_den = 14

# The singular values are s1 = r_y and s2 = r_x.
# The reciprocals are 1/s1 = 14/9 and 1/s2 = 14.

# Step 2: Use the formula for the Minkowski-Bouligand dimension for this class of fractals.
# The general formula is D = 1 + (log(N) - log(1/r_y)) / log(1/r_x)
# This can be simplified algebraically:
# D = 1 + (log(5) - log(14/9)) / log(14)
# D = (log(14) + log(5) - log(14) + log(9)) / log(14)
# D = (log(5) + log(9)) / log(14)
# D = log(5 * 9) / log(14)
# D = log(45) / log(14)

# Step 3: Calculate the numerical values for the final equation.
numerator_val = math.log(N * r_y_num) # log(5*9) = log(45)
denominator_val = math.log(r_x_den)  # log(14)
dimension = numerator_val / denominator_val

# Step 4: Print the result, showing the equation with the final numbers.
print("The Minkowski-Bouligand dimension D can be found using the simplified formula:")
print(f"D = (log({N}) + log({r_y_num})) / log({r_x_den})")
print(f"D = log({N * r_y_num}) / log({r_x_den})")
print(f"D = {numerator_val} / {denominator_val}")
print(f"D = {dimension}")

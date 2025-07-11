import math

# Step 1: Define the number of self-similar copies (N).
# This is the number of black keys in one octave.
N = 5

# Step 2: Define the scaling factor (r).
# The original keyboard is 3 units wide by 1 unit high.
# There are 7 white keys, so a white key's width is 3/7.
# A black key is half as wide, so its width is (1/2) * (3/7) = 3/14.
# A black key's height is 9/14 of the total height, so 9/14.
# To fit a 3r x r keyboard into a 3/14 x 9/14 space, the scaling
# factor 'r' is limited by the width (r <= 1/14) and height (r <= 9/14).
# We take the more restrictive factor.
r = 1 / 14

# The reciprocal of the scaling factor.
reciprocal_r = 1 / r

# Step 3: Calculate the Minkowski–Bouligand dimension using the formula D = log(N) / log(1/r).
dimension = math.log(N) / math.log(reciprocal_r)

# Step 4: Print the final equation and the result.
print("The formula for the dimension D is log(N) / log(1/r)")
print(f"The number of copies is N = {N}")
print(f"The scaling factor is r = 1/14, so 1/r = {int(reciprocal_r)}")
print(f"D = log({N}) / log({int(reciprocal_r)})")
print(f"The calculated dimension is: {dimension}")
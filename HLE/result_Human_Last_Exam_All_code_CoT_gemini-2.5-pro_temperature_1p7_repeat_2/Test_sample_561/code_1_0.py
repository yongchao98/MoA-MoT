import math

# The problem defines a self-affine fractal.
# N is the number of self-similar pieces at each iteration.
# In this case, it's the 5 black keys.
N = 5

# The overall keyboard is 3 units wide by 1 unit high.
# The transformation scales this overall size down to the size of a single black key.
# The scaling in the x-direction is sx = (black_key_width / octave_width) = ((3/7)/2) / 3 = 1/14
# The scaling in the y-direction is sy = (black_key_height / octave_height) = (9/14) / 1 = 9/14
sx = 1/14
sy = 9/14

# The box-counting dimension D can be found by evaluating the limit D = log(N(e))/log(1/e).
# For a self-affine fractal like this, the formula simplifies.
# Using a box size related to the smaller scaling factor sx, we find that the dimension is:
# D = log(N * (sy/sx)) / log(1/sx)
# where sy/sx = (9/14)/(1/14) = 9
# and 1/sx = 14
# So, N * (sy/sx) = 5 * 9 = 45.
# The formula becomes D = log(45) / log(14).

num_numerator = 45
num_denominator = 14

# Calculate the dimension
dimension = math.log(num_numerator) / math.log(num_denominator)

# Output the equation with each number and the final result.
print(f"The Minkowski-Bouligand dimension (D) is calculated by the equation:")
print(f"D = log({num_numerator}) / log({num_denominator})")
print(f"D = {dimension}")

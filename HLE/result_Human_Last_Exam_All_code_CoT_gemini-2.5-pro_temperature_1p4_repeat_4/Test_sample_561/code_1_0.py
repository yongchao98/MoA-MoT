import math

# Number of self-similar pieces (the number of black keys in an octave)
N = 5

# The horizontal scaling factor is the width of a black key (3/14) divided by the width of the keyboard (3).
# s_x = (3/14) / 3 = 1/14. The denominator of the scaling factor is 14.
s_x_denominator = 14

# The formula for the Minkowski-Bouligand dimension (D) of this self-affine system simplifies to:
# D = log(N) / log(1/s_x)
# which is log(N) / log(s_x_denominator)

# Calculate the dimension
dimension = math.log(N) / math.log(s_x_denominator)

# Print the equation and the result
print(f"The dimension D is calculated by the formula:")
print(f"D = log(Number of Copies) / log(1 / Horizontal Scaling Factor)")
print(f"D = log({N}) / log({s_x_denominator})")
print(f"D = {dimension}")
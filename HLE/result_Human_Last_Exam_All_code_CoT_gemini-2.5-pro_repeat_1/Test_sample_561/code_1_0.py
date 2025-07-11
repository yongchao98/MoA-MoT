import math

# This script calculates the Minkowski–Bouligand dimension for the described fractal.

# 1. Number of self-similar copies (N).
# A piano octave has 5 black keys. Each is replaced by a smaller keyboard.
N = 5

# 2. The scaling factor (r) is determined by fitting the original keyboard (3x1)
#    onto a black key ((3/14)x(9/14)). The scaling factor is r = 1/14.
#    The formula uses 1/r.
one_over_r = 14

# 3. Calculate the dimension using the formula D = log(N) / log(1/r).
dimension = math.log(N) / math.log(one_over_r)

# 4. Print the explanation, the final equation with its components, and the result.
print("The Minkowski–Bouligand dimension (D) is calculated from the number of self-similar copies (N) and the scaling factor (r).")
print(f"Number of copies (black keys), N = {N}")
print(f"The inverse of the scaling factor, 1/r = {one_over_r}")
print("\nThe final equation and result are:")
print(f"D = log({N}) / log({one_over_r}) = {dimension}")
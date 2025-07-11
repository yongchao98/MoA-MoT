import math

# Plan:
# 1. Define N, the number of self-similar copies (5 black keys).
# 2. Define s, the reciprocal of the scaling factor r. As derived in the plan,
#    the scaling is limited by the width, making r = 1/14, so s = 1/r = 14.
# 3. Apply the Minkowski–Bouligand dimension formula: D = log(N) / log(s)
# 4. Print the explanation, the final equation with values, and the result.

# Number of self-similar copies (number of black keys in an octave)
N = 5

# The reciprocal of the scaling factor (1/r)
# The scaling is constrained by fitting the 3-unit width of the keyboard
# into the 3/14-unit width of a black key. So, r = (3/14)/3 = 1/14.
# Thus, 1/r = 14.
s = 14

# Calculate the dimension
dimension = math.log(N) / math.log(s)

print("To find the Minkowski–Bouligand dimension, we use the formula: D = log(N) / log(1/r)")
print("-" * 20)
print(f"1. N, the number of self-similar copies, is the number of black keys: {N}")
print(f"2. 1/r, the reciprocal of the scaling factor, is determined by the geometry: {s}")
print("-" * 20)
print("Plugging the numbers into the formula, we get the final equation:")
print(f"D = log({N}) / log({s})")
print(f"\nThe calculated dimension is: {dimension}")
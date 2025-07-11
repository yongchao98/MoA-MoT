import math

# The number of self-similar copies (N) is the number of black keys in an octave.
N = 5

# The reciprocal of the scaling factor (1/r) is 14, as calculated in the plan.
one_over_r = 14

# The Minkowski–Bouligand dimension is D = log(N) / log(1/r)
dimension = math.log(N) / math.log(one_over_r)

print(f"The fractal is made of N = {N} self-similar copies.")
print(f"The reciprocal of the scaling factor is 1/r = {one_over_r}.")
print(f"The dimension D is calculated using the formula: D = log({N}) / log({one_over_r})")
print(f"The calculated Minkowski–Bouligand dimension is: {dimension}")

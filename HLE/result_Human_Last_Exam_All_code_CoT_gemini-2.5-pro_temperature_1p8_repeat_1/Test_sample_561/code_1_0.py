import math

# The number of self-similar pieces (N) is the number of black keys in an octave.
N = 5

# The inverse scaling factor (1/r) is determined by how much smaller the
# new keyboard is compared to the original. As calculated in the plan, r = 1/14.
one_over_r = 14

# The Minkowski–Bouligand dimension (D) for a self-similar set is calculated
# using the formula: D = log(N) / log(1/r)
dimension = math.log(N) / math.log(one_over_r)

print("Calculating the Minkowski–Bouligand dimension using the formula D = log(N) / log(1/r)")
print("-" * 70)
print(f"Number of self-similar copies (N): {N}")
print(f"Inverse scaling factor (1/r): {one_over_r}")
print("-" * 70)
print("The final equation is:")
print(f"D = log({N}) / log({one_over_r})")
print(f"D = {dimension}")
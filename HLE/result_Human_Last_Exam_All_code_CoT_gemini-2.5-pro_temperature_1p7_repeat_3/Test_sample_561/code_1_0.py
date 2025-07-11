import math

# Define the parameters based on the problem description.
# N is the number of self-similar copies, which is the number of black keys in an octave.
N = 5

# The larger scaling factor, s1, corresponds to the height reduction.
s1_num = 9
s1_den = 14
s1 = s1_num / s1_den

# The smaller scaling factor, s2, corresponds to the width reduction.
s2_num = 1
s2_den = 14
s2 = s2_num / s2_den

# The Minkowski-Bouligand dimension (D) for a self-affine fractal like this one
# is the solution to the equation: N * s1 * s2^(D - 1) = 1.

# We will now print out the equation with the specific numbers.
print("The dimension D is the solution to the following equation:")
print(f"{N} * ({s1_num}/{s1_den}) * ({s2_num}/{s2_den})^(D - 1) = 1")
print("")

# Solving this equation for D leads to: D = 1 + log(N * s1) / log(1/s2)
# D = 1 + log(5 * 9/14) / log(1 / (1/14))
# D = 1 + log(45/14) / log(14)

log_numerator_arg_num = N * s1_num
log_numerator_arg_den = s1_den
log_denominator_arg = s2_den

# Perform the final calculation.
dimension = 1 + math.log(log_numerator_arg_num / log_numerator_arg_den) / math.log(log_denominator_arg)

print("Solving for D, we get the expression:")
print(f"D = 1 + log({log_numerator_arg_num}/{log_numerator_arg_den}) / log({log_denominator_arg})")
print("")
print("The final calculated dimension is:")
print(dimension)
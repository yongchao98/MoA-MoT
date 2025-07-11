import math

# The problem is to find an upper bound for the maximum norm (k_k,inf) 
# in relation to the covolume (V) for quadratic number fields.
# This bound is derived from Minkowski's theorem.

# The constant in the inequality k_k_inf <= C * V depends on whether
# the field is real or imaginary.

# For real quadratic fields, the constant is 1/2.
constant_real = 1 / 2

# For imaginary quadratic fields, the constant is 4/pi.
constant_imaginary = 4 / math.pi

# An upper bound that is universally true for any quadratic field
# must use the larger of these two constants.
max_constant = max(constant_real, constant_imaginary)

# The larger constant is 4/pi, so this defines our upper bound.
# The final inequality is k_k_inf <= (4/pi) * V.

# As requested, we will output the numbers in the final equation.
numerator_in_bound = 4
denominator_in_bound_symbol = "pi"

print("The upper bound is derived from the larger of two cases (real and imaginary quadratic fields).")
print(f"The resulting inequality relating the maximum norm (k_k,inf) to the covolume (V) is:")
print(f"k_k_inf <= ( {numerator_in_bound} / {denominator_in_bound_symbol} ) * V")

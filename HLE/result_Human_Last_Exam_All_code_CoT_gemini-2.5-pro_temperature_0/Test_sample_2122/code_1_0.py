import math

# Based on the derivation, the temperature T is given by the formula:
# T(u) = -u + 2u / (1 + exp(2 - 2u^2))
# where u = cosh(2t).

# The target time is t_f = arccosh(2)/2.
# At this time, the value of u is:
# u_f = cosh(2 * arccosh(2)/2) = cosh(arccosh(2)) = 2.
u_f = 2

# We substitute u_f = 2 into the formula for T.
# The numbers in the final equation are:
val_u = 2
num_const = 2
exp_c1 = 2
exp_c2 = -2
exp_c3 = 2

print("The final equation for the temperature T is:")
print(f"T = -{val_u} + ({num_const} * {val_u}) / (1 + exp({exp_c1} + ({exp_c2}) * {val_u}^{exp_c3}))")

# Calculate the intermediate steps
exponent = exp_c1 + exp_c2 * (val_u ** exp_c3)
numerator = num_const * val_u
denominator = 1 + math.exp(exponent)
temperature = -val_u + numerator / denominator

print("\nSubstituting the value u = 2:")
print(f"T = -2 + 4 / (1 + exp(2 - 2 * 2^2))")
print(f"T = -2 + 4 / (1 + exp({exponent}))")
print(f"\nThe temperature at the specified time is: {temperature}")
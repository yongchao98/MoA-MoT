import math

# The problem has been analyzed and a parametric solution for the trajectory has been found.
# The position x is given by the function x(p), where p is the slope dy/dx.
# x(p) = -3 * (p + 1) / sqrt(2*p + 3)

# To find the position x0 where y = -3, we first need to find the corresponding
# value of the parameter p. This leads to the equation p^4 - 18p - 27 = 0.
# The physically relevant solution for p is 3.

p = 3

# We can now substitute p = 3 into the equation for x to find x0.

# Define the constants from the equation x(p)
coeff = -3
p_add = 1
p_mult = 2
const_add = 3

# Calculate x0
x0 = coeff * (p + p_add) / math.sqrt(p_mult * p + const_add)

# Output the equation with the numbers substituted
print("The equation for x0 with the value of p=3 is:")
print(f"x0 = {coeff} * ({p} + {p_add}) / sqrt({p_mult} * {p} + {const_add})")
print("\nStep-by-step calculation:")
print(f"x0 = {coeff * (p + p_add)} / sqrt({p_mult * p + const_add})")
print(f"x0 = {coeff * (p + p_add)} / sqrt({p_mult * p + const_add})")
print(f"x0 = {coeff * (p + p_add)} / {math.sqrt(p_mult * p + const_add)}")
print(f"x0 = {x0}")

# The final answer
print(f"\nThe position x0 where the particle reaches y = -3 is {x0}.")
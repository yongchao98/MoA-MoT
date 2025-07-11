import math

# Step 1: Define the value of the parameter p that corresponds to y = -3.
# This value was found by solving the equation p**4 - 18*p - 27 = 0.
p = 3

# Step 2: Define the terms of the equation for x0.
# The formula is x0 = -3 * (p + 1) / sqrt(2*p + 3).
coeff = -3
term1_numerator = p + 1
term1_denominator_inside_sqrt = 2 * p + 3

# Step 3: Calculate the numerator and the denominator.
numerator = coeff * term1_numerator
denominator = math.sqrt(term1_denominator_inside_sqrt)

# Step 4: Calculate the final position x0.
x0 = numerator / denominator

# Step 5: Print the calculation step-by-step.
print(f"The parameter p corresponding to the vertical coordinate y = -3 is p = {p}.")
print("The position x0 is found by substituting this value into the parametric equation for x:")
print(f"x0 = {coeff} * ({p} + 1) / sqrt(2*{p} + {3})")
print(f"x0 = {coeff * (p+1)} / sqrt({2*p+3})")
print(f"x0 = {numerator} / {denominator}")
print(f"Thus, the position x0 is {x0}.")

# This script calculates the maximum number of roots for t=5 based on the derived formula.

# Value of t for part (b)
t = 5

# The formula for the maximum number of roots is t * (t - 1) / 2.
# Let's calculate this for t = 5.

# Step 1: Calculate t - 1
t_minus_1 = t - 1

# Step 2: Calculate the numerator t * (t - 1)
numerator = t * t_minus_1

# Step 3: Calculate the final result by dividing by 2
# Using integer division as the result must be an integer.
result = numerator // 2

# Print the final equation with all the numbers, as requested.
print(f"The calculation for t = {t} is:")
print(f"({t} * ({t} - 1)) / 2 = ({t} * {t_minus_1}) / 2 = {numerator} / 2 = {result}")
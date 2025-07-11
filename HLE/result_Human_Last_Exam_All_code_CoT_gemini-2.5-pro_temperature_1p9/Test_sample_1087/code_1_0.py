import math

# Step 1: Define the constants involved in the formula.
val_5 = 5.0
val_2 = 2.0

# Step 2: Calculate the value of sqrt(5).
sqrt_5 = math.sqrt(val_5)

# Step 3: Calculate the numerator of the fraction inside the main square root.
numerator = val_5 - sqrt_5

# Step 4: Calculate the fraction inside the main square root.
fraction = numerator / val_2

# Step 5: Calculate the final value of r by taking the square root.
r = math.sqrt(fraction)

# Step 6: Print the formula and the final answer.
# The final formula for r is sqrt((5 - sqrt(5))/2)
print(f"The calculation is based on the formula: r = sqrt(({val_5} - sqrt({val_5})) / {val_2})")
print(f"The largest real number r is: {r}")

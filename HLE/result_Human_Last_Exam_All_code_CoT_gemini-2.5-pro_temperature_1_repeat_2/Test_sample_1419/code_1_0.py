import math

# This script derives and displays the leading order expression for the 
# Wilson-Fisher fixed point coupling u* in φ^4 theory near four dimensions.

# Step 1: Define the one-loop beta function in d = 4 - ϵ dimensions.
print("Step 1: The one-loop beta function for the coupling u is given by:")
print("β(u) = -ϵ*u + (3 * u^2) / (16 * π^2)\n")

# Step 2: Set the beta function to zero to find the fixed point u*.
print("Step 2: A non-trivial fixed point u* is found by setting β(u*) = 0:")
print("-ϵ*u* + (3 * (u*)^2) / (16 * π^2) = 0\n")

# Step 3: Solve the equation for u* (assuming u* ≠ 0).
print("Step 3: Solving for u* in terms of ϵ:")
print("Divide by u* -> -ϵ + (3 * u*) / (16 * π^2) = 0")
print("Rearrange -> (3 * u*) / (16 * π^2) = ϵ\n")

# Step 4: Present the final expression.
print("Step 4: The final leading order expression for the fixed point coupling is:")
numerator_coeff = 16
denominator_coeff = 3
pi_squared_str = "π^2"
epsilon_str = "ϵ"

# Output the final equation, showing each component clearly.
print(f"u* = ({numerator_coeff} * {pi_squared_str} / {denominator_coeff}) * {epsilon_str}")

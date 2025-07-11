# The Carman-Kozeny source term in the momentum equation for a two-phase
# melting simulation using the enthalpy-porosity method is often written as:
# S = -C * ((1 - f_l)^2 / (f_l^3 + b)) * u
# where 'C' is a large mushy zone constant, also referred to as the computational factor.

# In the prior published simulation-only work by Brent, Voller, and Reid (1988),
# before it was benchmarked against the melting of gallium, a specific value for C was used.
# This script identifies and prints that value.

# Define the components of the value from the paper.
base_value = 1.6
exponent = 6

# Calculate the final value.
computational_factor = base_value * (10**exponent)

# Print the final equation with each number.
# The format is chosen to match the answer choices.
print(f"The original computational factor (C) from the prior work was:")
print(f"{base_value} * 10^{exponent} = {computational_factor:.1e}")

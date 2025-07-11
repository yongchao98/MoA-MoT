# Define the given values for the firing level and membership degree
phi_k = 0.7
mu_g = 0.9

# A t-norm is used to calculate the rule activation level.
# We will use the minimum t-norm, which is the most common choice.
# The calculation is activation_level = min(phi_k, mu_g)
activation_level = min(phi_k, mu_g)

# Print the final equation and the result
print(f"The rule activation level using the minimum t-norm is calculated as follows:")
print(f"min({phi_k}, {mu_g}) = {activation_level}")
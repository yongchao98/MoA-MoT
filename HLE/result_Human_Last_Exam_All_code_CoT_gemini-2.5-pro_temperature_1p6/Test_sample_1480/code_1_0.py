# Define the given values
antecedent_firing_level = 0.7
consequent_membership_degree = 0.9

# A t-norm is used to combine these values. The most common t-norm is the minimum operator.
# The rule activation level is the result of this operation.
rule_activation_level = min(antecedent_firing_level, consequent_membership_degree)

# Print the final equation with the numbers and the result
print("The rule activation level is calculated as the t-norm (minimum) of the antecedent firing strength and the consequent membership degree.")
print(f"Rule Activation Level = min({antecedent_firing_level}, {consequent_membership_degree}) = {rule_activation_level}")
import numpy as np

# Given values for the fuzzy rule
# Antecedent activation level
antecedent_activation = 0.7

# Consequent membership degree
consequent_membership = 0.9

# A t-norm is used to combine these values. The most common t-norm in fuzzy
# logic is the minimum operator. We will use it to calculate the rule
# activation level.
# Rule Activation Level = T(antecedent_activation, consequent_membership)
# Using minimum t-norm: min(antecedent_activation, consequent_membership)
rule_activation_level = np.minimum(antecedent_activation, consequent_membership)

# Print the final equation with the values plugged in
print(f"The rule activation level is calculated by applying a t-norm (minimum operator) to the antecedent activation and the consequent membership.")
print(f"Result = min({antecedent_activation}, {consequent_membership})")
print(f"The final calculated rule activation level is: {rule_activation_level}")
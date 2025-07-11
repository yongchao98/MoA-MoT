# The given membership degrees for the rule's antecedents.
antecedent_1_membership = 0.7
antecedent_2_membership = 0.9

# The rule activation level is calculated by applying a t-norm to the
# antecedent membership degrees. The most common t-norm is the minimum operator.
activation_level = min(antecedent_1_membership, antecedent_2_membership)

# Print the final result showing the full equation.
print(f"The rule activation level using the minimum t-norm is calculated as follows:")
print(f"min({antecedent_1_membership}, {antecedent_2_membership}) = {activation_level}")
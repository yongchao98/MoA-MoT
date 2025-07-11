# The given membership values for the antecedents of the fuzzy rule.
# Phi_k(x') is interpreted as the membership degree of the first antecedent.
antecedent_1_membership = 0.7
# μG_k_j(y_j) is interpreted as the membership degree of the second antecedent.
antecedent_2_membership = 0.9

# A t-norm operation is used to combine the antecedent memberships to find the
# overall rule activation level. We will use the minimum t-norm, which is the
# most common choice.
rule_activation_level = min(antecedent_1_membership, antecedent_2_membership)

# Display the final calculation and result.
print(f"The rule activation level is calculated by applying a t-norm to the antecedent membership degrees.")
print(f"Using the minimum t-norm, the equation is:")
print(f"min({antecedent_1_membership}, {antecedent_2_membership}) = {rule_activation_level}")
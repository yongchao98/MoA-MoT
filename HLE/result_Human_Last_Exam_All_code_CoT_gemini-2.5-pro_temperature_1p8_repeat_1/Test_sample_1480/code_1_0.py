# Define the membership values from the rule's antecedents.
# Phi_k(x') can be interpreted as the membership value for the first condition.
val1 = 0.7
# muG_k_j(y_j) can be interpreted as the membership value for the second condition.
val2 = 0.9

# A t-norm is used to combine the membership values for a fuzzy 'AND' operation.
# The most common t-norm is the minimum operator.
# We calculate the rule activation level by taking the minimum of the two values.
activation_level = min(val1, val2)

# Print the final equation and the result.
print(f"Using the Minimum t-norm, the rule activation level is: min({val1}, {val2}) = {activation_level}")
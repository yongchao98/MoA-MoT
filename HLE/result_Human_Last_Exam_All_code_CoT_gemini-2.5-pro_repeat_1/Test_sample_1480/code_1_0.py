import numpy as np

# Step 1: Define the given membership values from the antecedent of the rule.
# We interpret the provided values, 0.7 and 0.9, as the membership grades to be combined.
mu1 = 0.7
mu2 = 0.9

# Step 2: Apply the t-norm operation.
# The most common t-norm is the minimum operator, which we will use here.
# This represents the fuzzy 'AND' operation.
activation_level = np.minimum(mu1, mu2)

# Step 3: Print the final equation and the result.
# The output shows the operation performed and its result.
print(f"The rule activation level is the t-norm of the antecedent membership values.")
print(f"Using the minimum t-norm: min({mu1}, {mu2}) = {activation_level}")

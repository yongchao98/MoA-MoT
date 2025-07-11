# Define the given values
antecedent_firing_level = 0.7
consequent_membership = 0.9

# A t-norm is a fuzzy logic 'AND' operator. The most common t-norm is the minimum function.
# Although the problem mentions an Interval Type-3 Fuzzy Set, the provided values are scalars,
# so the calculation is equivalent to a Type-1 fuzzy logic operation.
# We apply the minimum t-norm to the antecedent firing level and the consequent membership degree.
activation_level = min(antecedent_firing_level, consequent_membership)

# Print the explanation and the result
print("The rule activation level is calculated using the minimum t-norm, which finds the minimum of the antecedent firing level and the consequent membership degree.")
print(f"Calculation: min({antecedent_firing_level}, {consequent_membership}) = {activation_level}")
# The given values for the membership grades in the rule's antecedent.
# While the problem mentions an Interval Type-3 Fuzzy Set, the inputs provided are
# crisp numbers, simplifying the calculation to a standard t-norm operation.
value1 = 0.7
value2 = 0.9

# A t-norm is a fuzzy logic operator for 'AND'. The most common t-norm is the minimum function.
# We will use the minimum t-norm to calculate the rule activation level.
activation_level = min(value1, value2)

print("To find the rule activation level, we combine the membership values of the antecedent using a t-norm.")
print("Using the minimum t-norm (the most common choice):")
print(f"Activation Level = min({value1}, {value2})")
print(f"The resulting rule activation level is: {activation_level}")
import math

# Define the parameters from the problem description for ontology O1
b = 3  # Branching factor
d = 4  # Depth

# According to the correct theoretical model (Statement D), the completeness
# is bounded by the ratio of the natural logarithm of the branching factor
# to the depth of the hierarchy.
# This formula models the dilution of structural information over the hierarchy's depth.
completeness_bound = math.log(b) / d

# Print the final equation with the numbers plugged in, as requested.
# The result shows the calculated upper bound on theoretical completeness.
print(f"The theoretical completeness of the mapping is bounded by ln(b)/d.")
print(f"Using the given values:")
print(f"ln({b}) / {d} = {completeness_bound:.4f}")
print("\nThis value represents the theoretical upper limit on how much of O₁'s logical structure can be preserved in O₂.")

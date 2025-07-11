import math

# Based on the analysis of the polynomial's coefficients using Vieta's formulas,
# the four roots of the equation were determined to be:
# sqrt(14), 2*sqrt(6), sqrt(34), and 2*sqrt(11).

# Define the four roots.
roots = [
    math.sqrt(14),
    2 * math.sqrt(6),
    math.sqrt(34),
    2 * math.sqrt(11)
]

# Sort the roots in increasing order.
roots.sort()

# Print the sorted roots with clear labeling.
print("The four roots of the polynomial in increasing order are:")
for i, root in enumerate(roots):
    print(f"x_{i+1} = {root}")

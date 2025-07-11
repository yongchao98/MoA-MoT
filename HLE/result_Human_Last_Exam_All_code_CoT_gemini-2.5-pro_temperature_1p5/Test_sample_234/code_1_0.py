# This script formalizes the counting of properties that are always true.

# A dictionary mapping each property to a boolean indicating if it must always be true.
property_analysis = {
    "Open": True,
    "Closed": False,
    "Connected": False,
    "Compact": False,
    "Dense": True,
    "Connected complement": False,
    "Trivial first singular homology group": False
}

# Calculate the total count of properties that are always true.
# The sum() function on a list of booleans treats True as 1 and False as 0.
always_true_count = sum(property_analysis.values())

# Create a representation of the summation for the output, as requested.
# It will be a string like "1 + 1" if the count is 2.
summands = ["1"] * always_true_count
equation = " + ".join(summands)

# Print the final result and the explanatory "equation".
print("Based on the mathematical analysis of the seven properties:")
print(f"Number of properties that must always be true: {always_true_count}")
if always_true_count > 0:
    print(f"This is the sum of the properties found to be true: {equation} = {always_true_count}")
else:
    # This case would apply if no properties were always true.
    print("This is the sum of the properties found to be true: 0")

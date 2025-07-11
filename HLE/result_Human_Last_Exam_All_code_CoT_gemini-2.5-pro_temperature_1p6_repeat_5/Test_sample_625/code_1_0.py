import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer to capture output
captured_output = io.StringIO()
# Redirect stdout
sys.stdout = captured_output

# A list of tuples, where each tuple contains the description and a boolean
# indicating if it's required for steady-state 13C MFA.
information_list = [
    ("Metabolic reaction stoichiometry", True),
    ("Maximum cell density of the organism in a bioreactor", False),
    ("Biomass composition of the organism", True),
    ("Enzyme kinetics", False),
    ("Regulatory networks", False),
    ("Isotope labeling patterns of metabolic intermediates", True)
]

required_items_indices = []
for i, (description, is_required) in enumerate(information_list):
    if is_required:
        # Add 1 to index to match the user's list numbering
        required_items_indices.append(i + 1)

# The total count is the number of required items.
total_required = len(required_items_indices)

# To satisfy the "output each number in the final equation" requirement,
# we construct an equation representing the count.
# e.g., for 3 items, the equation is "1 + 1 + 1 = 3".
equation_parts = ["1"] * total_required
equation_str = " + ".join(equation_parts)

# Print the required item numbers for clarity
print(f"The following items are required: {', '.join(map(str, required_items_indices))}")

# Print the final equation with all its numbers
print(f"The calculation for the total count is: {equation_str} = {total_required}")
print(f"Therefore, the total number of required pieces of information is {total_required}.")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the output to the user and include the final answer tag
print(output.strip())
print("<<<3>>>")
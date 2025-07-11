import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# This script calculates the number of required pieces of information for 13C MFA.

# A list of boolean flags representing the requirements for the 6 items.
# True means required, False means not required.
requirements_flags = [
    True,   # 1. Metabolic reaction stoichiometry
    False,  # 2. Maximum cell density
    True,   # 3. Biomass composition of the organism
    False,  # 4. Enzyme kinetics
    False,  # 5. Regulatory networks
    True    # 6. Isotope labeling patterns of metabolic intermediates
]

# A list of the descriptions for context in the output
requirements_descriptions = [
    "Metabolic reaction stoichiometry",
    "Maximum cell density of the organism in a bioreactor",
    "Biomass composition of the organism (i.e. fraction of protein, lipids, and carbohydrates)",
    "Enzyme kinetics",
    "Regulatory networks",
    "Isotope labeling patterns of metabolic intermediates"
]


# Identify which items are required
required_items_indices = [i + 1 for i, req in enumerate(requirements_flags) if req]
count_of_required_items = len(required_items_indices)

# Build an "equation" string to show the calculation
# e.g., "1 + 1 + 1" for three required items
equation_parts = ["1"] * count_of_required_items
equation_str = " + ".join(equation_parts)


print("The following information is required for 13C MFA at steady state:")
for i in required_items_indices:
    print(f"{i}. {requirements_descriptions[i-1]}")

print("\nCalculating the total count:")
# Print the final calculation
print(f"{equation_str} = {count_of_required_items}")

# Restore original stdout
sys.stdout = original_stdout
# Get the content from the buffer
output = captured_output.getvalue()

# Print the captured output to the user
print(output)

# Final answer in the specified format
final_answer = count_of_required_items
print(f"<<<{final_answer}>>>")
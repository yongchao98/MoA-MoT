# List of information items and whether they are required for 13C MFA.
# The tuple format is (description, is_required_boolean).
information_list = [
    ("Metabolic reaction stoichiometry", True),
    ("Maximum cell density of the organism in a bioreactor", False),
    ("Biomass composition of the organism", True),
    ("Enzyme kinetics", False),
    ("Regulatory networks", False),
    ("Isotope labeling patterns of metabolic intermediates", True)
]

# Initialize a counter for the required items.
required_item_count = 0
required_items_names = []

print("Analyzing the requirements for 13C Metabolic Flux Analysis...\n")

# Iterate through the list to check each item.
for item, is_required in information_list:
    if is_required:
        status = "Required"
        required_item_count += 1
        required_items_names.append(item)
    else:
        status = "Not Required"
    print(f"- {item}: {status}")

print("\n-------------------------------------------------")
print("The required pieces of information are:")
for name in required_items_names:
    print(f"  - {name}")

print("\nFinal Count:")
# The instruction "output each number in the final equation" is interpreted
# as showing the final calculated number clearly.
print(f"Total number of required items = {required_item_count}")
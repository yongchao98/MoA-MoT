import collections

# Define the partial charges for methanol based on the selected reasonable model.
# The chosen model is charge-neutral and respects chemical principles.
# Atom types: C (Carbon), Hc (Methyl Hydrogen), O (Oxygen), Ho (Hydroxyl Hydrogen)
charges = {
    "Carbon": 0.1450,
    "Methyl H1": 0.0400,
    "Methyl H2": 0.0400,
    "Methyl H3": 0.0400,
    "Oxygen": -0.6830,
    "Hydroxyl H": 0.4180,
}

# The atoms in a methanol molecule
atoms = ["Carbon", "Methyl H1", "Methyl H2", "Methyl H3", "Oxygen", "Hydroxyl H"]

# Print the proposed charge for each atom
print("Proposed Partial Charges for Methanol (CH3OH):")
for atom_name, charge_value in charges.items():
    print(f"{atom_name}:\t{charge_value:.4f}")

# Verify that the total charge of the molecule is zero
total_charge = sum(charges.values())

print("\nVerifying a total charge of zero:")
# Build the equation string dynamically
equation_parts = []
for atom_name in atoms:
    charge_value = charges[atom_name]
    # Format the number to show its sign
    equation_parts.append(f"{charge_value:+.4f}")

# Join parts with ' + ' but handle the leading plus sign of the next number
equation_str = " + ".join(equation_parts).replace("+ -", "- ")

print(f"{equation_str} = {total_charge:.4f}")
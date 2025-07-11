# Proposed partial charges for a new all-atom model of methanol (CH3OH).
charges = {
    'Carbon': 0.1450,
    'Methyl H1': 0.0400,
    'Methyl H2': 0.0400,
    'Methyl H3': 0.0400,
    'Oxygen': -0.6830,
    'Hydroxyl H': 0.4180
}

# Calculate the total charge to verify neutrality.
total_charge = sum(charges.values())

# Print the proposed charge for each atom.
print("Proposed partial charge assignments for Methanol:")
for atom, charge in charges.items():
    # The format ensures trailing zeros are printed, matching the input format.
    print(f"{atom}:\t{charge:.4f}")

# Print the final equation summing the charges
charge_values = list(charges.values())
equation_parts = [f"({c:.4f})" if c < 0 else f"{c:.4f}" for c in charge_values]
equation = " + ".join(equation_parts)
print(f"\nSum of charges: {equation} = {total_charge:.4f}")

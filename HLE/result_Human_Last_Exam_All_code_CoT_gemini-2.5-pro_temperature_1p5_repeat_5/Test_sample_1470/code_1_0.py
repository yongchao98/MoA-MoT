# Propose a set of reasonable partial charges for a methanol (CH3OH) molecule.
# The charges must sum to zero to ensure the molecule is neutral.
# We will check the validity of the charges from the most chemically reasonable option.

charges = {
    'Carbon': 0.1450,
    'Methyl H1': 0.0400,
    'Methyl H2': 0.0400,
    'Methyl H3': 0.0400,
    'Oxygen': -0.6830,
    'Hydroxyl hydrogen': 0.4180
}

# Print the proposed charges
print("Proposed Partial Charges for Methanol (in fundamental charge units):")
for atom, charge in charges.items():
    print(f"{atom}:\t{charge:.4f}")

# Verify that the sum of charges is zero
total_charge = sum(charges.values())

# Print the neutrality check equation
print("\nVerifying Molecular Neutrality:")
charge_values = list(charges.values())
equation_str = " + ".join([f"({c:.4f})" for c in charge_values])
print(f"Sum of charges: {equation_str} = {total_charge:.4f}")

if round(total_charge, 4) == 0.0:
    print("\nThe molecule is electrically neutral. This is a valid and reasonable set of charges.")
else:
    print("\nWarning: The molecule is not electrically neutral. This is not a valid set of charges.")

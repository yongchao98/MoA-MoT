import collections

# Define the proposed partial charges based on the most reasonable choice.
# These charges are chemically sound, respect symmetry, and sum to zero.
charges = collections.OrderedDict([
    ("Carbon", 0.1450),
    ("Methyl H1", 0.0400),
    ("Methyl H2", 0.0400),
    ("Methyl H3", 0.0400),
    ("Oxygen", -0.6830),
    ("Hydroxyl H", 0.4180)
])

# Calculate the total charge
total_charge = sum(charges.values())

# Print the proposed assignments
print("Proposed partial charge assignments for methanol:")
for atom, charge in charges.items():
    # The ' ' prefix aligns the positive numbers with the negative ones
    print(f"{atom}:\t{' ' if charge >= 0 else ''}{charge:.4f}")

print("\nVerification of net charge neutrality:")
# Build the equation string
charge_values = [
    charges["Carbon"], 
    charges["Methyl H1"], 
    charges["Methyl H2"], 
    charges["Methyl H3"], 
    charges["Oxygen"], 
    charges["Hydroxyl H"]
]
equation_parts = [f"({val:.4f})" if val < 0 else f"{val:.4f}" for val in charge_values]
equation_str = " + ".join(equation_parts)

print(f"Sum of charges: {equation_str} = {total_charge:.4f}")

# Final confirmation
if abs(total_charge) < 1e-9:
    print("\nThe sum of the partial charges is zero, confirming the molecule is neutral.")
else:
    print(f"\nWarning: The total charge is {total_charge:.4f}, not zero.")

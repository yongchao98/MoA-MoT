# Define the partial charges for the atoms in methanol
charges = {
    'Carbon': 0.1450,
    'Methyl H1': 0.0400,
    'Methyl H2': 0.0400,
    'Methyl H3': 0.0400,
    'Oxygen': -0.6830,
    'Hydroxyl hydrogen': 0.4180
}

# Print the proposed charges
print("Proposed Partial Charges for Methanol:")
for atom, charge in charges.items():
    print(f"{atom}:\t{charge:.4f}")

# Verify that the total charge is zero
total_charge = sum(charges.values())

# Create the equation string for verification
charge_values = [charges['Carbon'], charges['Methyl H1'], charges['Methyl H2'], charges['Methyl H3'], charges['Oxygen'], charges['Hydroxyl hydrogen']]
equation_str = " + ".join(map(str, charge_values))

print("\nVerification of Net Charge Neutrality:")
print(f"Sum = {equation_str} = {total_charge:.4f}")

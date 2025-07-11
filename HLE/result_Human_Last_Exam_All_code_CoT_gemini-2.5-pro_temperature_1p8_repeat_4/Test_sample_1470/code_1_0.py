import collections

# Proposing reasonable partial charges for a new all-atom model of methanol.
# Based on analysis, the parameters from choice D are the most reasonable.
proposed_charges = collections.OrderedDict([
    ('Carbon', 0.1450),
    ('Methyl H1', 0.0400),
    ('Methyl H2', 0.0400),
    ('Methyl H3', 0.0400),
    ('Oxygen', -0.6830),
    ('Hydroxyl hydrogen', 0.4180)
])

print("Proposed Partial Charges for Methanol (CH3OH):")
for atom, charge in proposed_charges.items():
    # Using '{: <18}' to align the charge values
    print(f"{atom + ':': <18} {charge: .4f}")

# Verify that the sum of the charges is zero (net neutrality).
print("\nVerifying Net Neutrality:")
total_charge = sum(proposed_charges.values())

# Create a string representation of the sum equation
# Explicitly including each number for the final equation.
charge_values = list(proposed_charges.values())
equation_string = f"({charge_values[0]:.4f}) + ({charge_values[1]:.4f}) + ({charge_values[2]:.4f}) + ({charge_values[3]:.4f}) + ({charge_values[4]:.4f}) + ({charge_values[5]:.4f})"

print("Sum of charges calculation:")
print(f"{equation_string} = {total_charge:.4f}")

if abs(total_charge) < 1e-9:
    print("\nThe proposed charges result in a neutral molecule, which is a requirement for the model.")
else:
    print(f"\nWarning: The proposed charges result in a non-neutral molecule with a net charge of {total_charge:.4f}.")

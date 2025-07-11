import collections

# Based on chemical principles, option D is the most reasonable set of charges.
# Let's define the charges for each atom in the methanol molecule.
# C: Carbon
# H_m1, H_m2, H_m3: Methyl Hydrogens
# O: Oxygen
# H_o: Hydroxyl Hydrogen
charges = collections.OrderedDict([
    ('Carbon', 0.1450),
    ('Methyl H1', 0.0400),
    ('Methyl H2', 0.0400),
    ('Methyl H3', 0.0400),
    ('Oxygen', -0.6830),
    ('Hydroxyl hydrogen', 0.4180)
])

# Print the proposed charges
print("Proposed Partial Charges for Methanol (CH3OH):")
for atom, charge in charges.items():
    # The format specifier ensures trailing zeros are printed
    print(f"{atom}:\t{charge:.4f}")

# Verify that the total charge is zero and print the calculation
charge_values = list(charges.values())
total_charge = sum(charge_values)

print("\nVerifying charge neutrality:")
# Build the equation string, handling negative numbers
equation_parts = [f"{c:.4f}" for c in charge_values]
equation_str = " + ".join(equation_parts).replace("+ -", "- ")
print(f"Sum = {equation_str}")
print(f"Total Charge = {total_charge:.4f}")

if round(total_charge, 4) == 0.0:
    print("\nThe molecule is charge-neutral, as required.")
else:
    print("\nWarning: The molecule is not charge-neutral.")

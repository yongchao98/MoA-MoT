import collections

# Based on chemical principles, Choice D provides a reasonable set of partial charges.
# 1. The charges reflect the electronegativity differences (O << C < H).
# 2. The three chemically equivalent methyl hydrogens have identical charges.
# 3. The net charge of the molecule sums to zero.

charges = collections.OrderedDict([
    ("Carbon", 0.1450),
    ("Methyl H1", 0.0400),
    ("Methyl H2", 0.0400),
    ("Methyl H3", 0.0400),
    ("Oxygen", -0.6830),
    ("Hydroxyl hydrogen", 0.4180)
])

print("Proposed Partial Charges for Methanol:")
for atom, charge in charges.items():
    # The format specifier ensures trailing zeros are printed
    print(f"{atom}:\t{charge:.4f}")

# Verify that the total charge is zero and print the equation
total_charge = sum(charges.values())
charge_values = list(charges.values())

print("\nVerification of Net Charge (Sum = 0):")
# We build the equation string, handling negative numbers correctly
equation_parts = [f"({c:.4f})" for c in charge_values]
equation_str = " + ".join(equation_parts)

print(f"{equation_str} = {total_charge:.4f}")

print("\n<<<D>>>")
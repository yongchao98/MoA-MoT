# Propose a reasonable set of partial charges for a new all-atom model of methanol.

# The selected charges are based on chemical principles:
# 1. The molecule must be neutral (sum of charges = 0).
# 2. Charges should reflect electronegativity (O is highly negative, H attached to O is highly positive).
# 3. Chemically equivalent atoms (the three methyl hydrogens) should have the same charge.
# Choice D fits these criteria best.

charges = {
    "Carbon": 0.1450,
    "Methyl H1": 0.0400,
    "Methyl H2": 0.0400,
    "Methyl H3": 0.0400,
    "Oxygen": -0.6830,
    "Hydroxyl H": 0.4180,
}

print("Proposed Partial Charges for Methanol (in fundamental charge units):\n")

total_charge = 0
for atom, charge in charges.items():
    print(f"{atom}:\t{charge:.4f}")
    total_charge += charge

print("\n-----------------------------")
print(f"Total molecular charge: {total_charge:.4f}")
print("-----------------------------")

print("\nFinal Equation:")
charge_values = list(charges.values())
equation_str = " + ".join([f"({c:.4f})" for c in charge_values])
print(f"{equation_str} = {total_charge:.4f}")

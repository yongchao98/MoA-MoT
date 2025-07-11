# Propose reasonable partial charges for a new all-atom model of methanol.

# Atom partial charges in fundamental charge units (e).
charge_C = 0.1450
charge_H_methyl = 0.0400
charge_O = -0.6830
charge_H_hydroxyl = 0.4180

# The three methyl hydrogens are chemically equivalent and have the same charge.
charge_H1 = charge_H_methyl
charge_H2 = charge_H_methyl
charge_H3 = charge_H_methyl

# Calculate the total charge to ensure the molecule is neutral.
total_charge = charge_C + charge_H1 + charge_H2 + charge_H3 + charge_O + charge_H_hydroxyl

print("Proposed Partial Charge Assignments for Methanol (CH3OH):")
print("-" * 55)
print(f"Carbon (C):\t\t{charge_C:.4f}")
print(f"Methyl Hydrogen (H1):\t{charge_H1:.4f}")
print(f"Methyl Hydrogen (H2):\t{charge_H2:.4f}")
print(f"Methyl Hydrogen (H3):\t{charge_H3:.4f}")
print(f"Oxygen (O):\t\t{charge_O:.4f}")
print(f"Hydroxyl Hydrogen (H_O):\t{charge_H_hydroxyl:.4f}")
print("-" * 55)

# Print the equation verifying the total charge is zero.
print("\nVerification of Net Charge:")
print(f"({charge_C:.4f}) + ({charge_H1:.4f}) + ({charge_H2:.4f}) + ({charge_H3:.4f}) + ({charge_O:.4f}) + ({charge_H_hydroxyl:.4f}) = {total_charge:.4f}")

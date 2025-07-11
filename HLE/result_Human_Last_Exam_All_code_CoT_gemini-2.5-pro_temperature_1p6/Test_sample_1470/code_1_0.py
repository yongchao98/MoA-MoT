# Propose and verify a reasonable set of partial charges for methanol (CH3OH).

# Atom charges from the selected answer choice.
# C(H3)-O-H
charge_C = 0.1450
charge_H_methyl = 0.0400 # Charge for each of the three methyl hydrogens
charge_O = -0.6830
charge_H_hydroxyl = 0.4180

# The three methyl hydrogens are chemically equivalent and should have the same charge.
charge_H1 = charge_H_methyl
charge_H2 = charge_H_methyl
charge_H3 = charge_H_methyl

# Calculate the total charge of the molecule to verify neutrality.
total_charge = charge_C + charge_H1 + charge_H2 + charge_H3 + charge_O + charge_H_hydroxyl

# Print the charge assignments.
print("Proposed Partial Charges for Methanol (CH3OH):")
print(f"Carbon: {charge_C:.4f}")
print(f"Methyl H1: {charge_H1:.4f}")
print(f"Methyl H2: {charge_H2:.4f}")
print(f"Methyl H3: {charge_H3:.4f}")
print(f"Oxygen: {charge_O:.4f}")
print(f"Hydroxyl Hydrogen: {charge_H_hydroxyl:.4f}")
print("\nVerifying the total charge of the molecule is zero:")
# The format specifier ensures all numbers are printed with consistent decimal places.
print(f"{charge_C:.4f} (C) + {charge_H1:.4f} (H) + {charge_H2:.4f} (H) + {charge_H3:.4f} (H) + {charge_O:.4f} (O) + {charge_H_hydroxyl:.4f} (H) = {total_charge:.4f}")
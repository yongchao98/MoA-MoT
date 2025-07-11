# Proposing partial charges for a new methanol (CH3OH) model.
# The chosen set must be chemically reasonable and sum to zero for a neutral molecule.

# Charges from Option D
charge_C = 0.1450
charge_H_methyl = 0.0400  # For each of the 3 methyl hydrogens
charge_O = -0.6830
charge_H_hydroxyl = 0.4180

# The methanol molecule has 1 Carbon, 3 methyl Hydrogens, 1 Oxygen, and 1 hydroxyl Hydrogen.
total_charge = charge_C + (3 * charge_H_methyl) + charge_O + charge_H_hydroxyl

print("Proposed Partial Charges for Methanol (CH3OH):")
print(f"Carbon: {charge_C:.4f}")
print(f"Methyl H1: {charge_H_methyl:.4f}")
print(f"Methyl H2: {charge_H_methyl:.4f}")
print(f"Methyl H3: {charge_H_methyl:.4f}")
print(f"Oxygen: {charge_O:.4f}")
print(f"Hydroxyl hydrogen: {charge_H_hydroxyl:.4f}")
print("\nVerifying Net Neutrality:")
# The problem asks to output each number in the final equation.
print(f"Sum = {charge_C:.4f} + (3 * {charge_H_methyl:.4f}) + ({charge_O:.4f}) + {charge_H_hydroxyl:.4f}")
print(f"Total Charge = {total_charge:.4f}")

if round(total_charge, 4) == 0.0:
    print("\nThe charges sum to zero. This is a valid set for a neutral molecule.")
else:
    print("\nWarning: The charges do not sum to zero.")

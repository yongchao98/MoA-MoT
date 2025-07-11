import sys

# This script proposes a reasonable set of partial charges for methanol (CH3OH)
# and verifies that the molecule is electrically neutral.

# Define the partial charges for each atom based on chemical principles.
# These values are chosen from the provided options as the most reasonable set.
charge_C = 0.1450
charge_H_methyl = 0.0400
charge_O = -0.6830
charge_H_hydroxyl = 0.4180

# For clarity, assign charges to each of the three equivalent methyl hydrogens.
charge_H1 = charge_H_methyl
charge_H2 = charge_H_methyl
charge_H3 = charge_H_methyl

# Calculate the total charge of the molecule.
total_charge = charge_C + charge_H1 + charge_H2 + charge_H3 + charge_O + charge_H_hydroxyl

# Print the proposed charge assignments for the user.
print("Proposed partial charges for a new methanol (CH3OH) model:")
print(f"{'Atom':<20} {'Charge (e)':<10}")
print("-" * 30)
print(f"{'Carbon':<20} {charge_C: 9.4f}")
print(f"{'Methyl H1':<20} {charge_H1: 9.4f}")
print(f"{'Methyl H2':<20} {charge_H2: 9.4f}")
print(f"{'Methyl H3':<20} {charge_H3: 9.4f}")
print(f"{'Oxygen':<20} {charge_O: 9.4f}")
print(f"{'Hydroxyl hydrogen':<20} {charge_H_hydroxyl: 9.4f}")
print("-" * 30)

# Print the verification equation to show the molecule is neutral.
# The 's' format specifier adds a '+' for positive numbers for clarity.
print("\nVerification of electroneutrality:")
print(f"({charge_C:+.4f}) + ({charge_H1:+.4f}) + ({charge_H2:+.4f}) + ({charge_H3:+.4f}) + ({charge_O:+.4f}) + ({charge_H_hydroxyl:+.4f}) = {total_charge:.4f}")

# The format function in python 3.8+ for f-strings can handle this more cleanly
# For compatibility, we'll stick to the above method.
# For example, in Py 3.8+:
# print(f"{charge_C=:.4f} + {charge_H1=:.4f} ...")

# Check if the total charge is close to zero within a small tolerance.
if abs(total_charge) < 1e-9:
    print("\nThe molecule is electrically neutral.")
else:
    print(f"\nWarning: The molecule is not neutral. Total charge is {total_charge:.4f}.")

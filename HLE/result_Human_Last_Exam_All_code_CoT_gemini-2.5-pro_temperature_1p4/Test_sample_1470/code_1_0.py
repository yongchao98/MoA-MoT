# Propose a reasonable set of partial charges for methanol (CH3OH).
# This code adopts the physically sound charge distribution from Choice D.

# Assign charges to each atom type.
# The three methyl hydrogens are equivalent and share the same charge.
charge_C = 0.1450
charge_H_methyl = 0.0400
charge_O = -0.6830
charge_H_hydroxyl = 0.4180

# The molecule consists of 1 Carbon, 3 methyl Hydrogens, 1 Oxygen, and 1 hydroxyl Hydrogen.
total_charge = charge_C + (3 * charge_H_methyl) + charge_O + charge_H_hydroxyl

# Print the proposed charge for each atom in the molecule.
print("Proposed Partial Charges for Methanol (in fundamental charge units):")
print(f"Carbon:            {charge_C:.4f}")
print(f"Methyl Hydrogen 1:   {charge_H_methyl:.4f}")
print(f"Methyl Hydrogen 2:   {charge_H_methyl:.4f}")
print(f"Methyl Hydrogen 3:   {charge_H_methyl:.4f}")
print(f"Oxygen:            {charge_O:.4f}")
print(f"Hydroxyl Hydrogen:   {charge_H_hydroxyl:.4f}")
print("-" * 35)

# Verify charge neutrality by printing the full summation equation.
# This fulfills the requirement to output each number in the final equation.
print("Verification of Net Charge Neutrality:")
print(f"Equation: {charge_C} + {charge_H_methyl} + {charge_H_methyl} + {charge_H_methyl} + ({charge_O}) + {charge_H_hydroxyl} = {total_charge:.4f}")
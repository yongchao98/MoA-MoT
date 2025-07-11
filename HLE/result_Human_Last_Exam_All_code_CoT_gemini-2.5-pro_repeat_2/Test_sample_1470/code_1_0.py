# Proposed partial charges for methanol (CH3OH) based on chemical principles.
# The molecule consists of 1 Carbon, 3 methyl Hydrogens, 1 Oxygen, and 1 hydroxyl Hydrogen.

charge_C = 0.1450
charge_H_methyl = 0.0400
charge_O = -0.6830
charge_H_hydroxyl = 0.4180

# The three methyl hydrogens are chemically equivalent and should have the same charge.
num_H_methyl = 3

# Calculate the total charge of the molecule
total_charge = charge_C + (num_H_methyl * charge_H_methyl) + charge_O + charge_H_hydroxyl

# Print the proposed charge assignments
print("Proposed Partial Charge Assignments for Methanol:")
print(f"Carbon: {charge_C:.4f}")
print(f"Methyl H1: {charge_H_methyl:.4f}")
print(f"Methyl H2: {charge_H_methyl:.4f}")
print(f"Methyl H3: {charge_H_methyl:.4f}")
print(f"Oxygen: {charge_O:.4f}")
print(f"Hydroxyl hydrogen: {charge_H_hydroxyl:.4f}")
print("\nVerifying Net Charge Neutrality:")
print(f"Equation: C + 3*H(methyl) + O + H(hydroxyl)")
print(f"Calculation: {charge_C:.4f} + 3*({charge_H_methyl:.4f}) + ({charge_O:.4f}) + {charge_H_hydroxyl:.4f} = {total_charge:.4f}")

if round(total_charge, 8) == 0:
    print("\nThe sum of the partial charges is 0. This is a valid charge set.")
else:
    print(f"\nWarning: The sum of the partial charges is {total_charge}, not 0. This is not a valid charge set.")

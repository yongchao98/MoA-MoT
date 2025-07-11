# Define the partial charges for each atom in the methanol molecule (CH3OH)
# based on a chemically reasonable model.
charge_C = 0.1450
charge_H_methyl = 0.0400
charge_O = -0.6830
charge_H_hydroxyl = 0.4180

# The methanol molecule has 1 Carbon, 3 methyl Hydrogens, 1 Oxygen, and 1 hydroxyl Hydrogen.
num_C = 1
num_H_methyl = 3
num_O = 1
num_H_hydroxyl = 1

# Calculate the total charge of the molecule
total_charge = (num_C * charge_C) + \
               (num_H_methyl * charge_H_methyl) + \
               (num_O * charge_O) + \
               (num_H_hydroxyl * charge_H_hydroxyl)

# Print the proposed partial charges for each atom
print("Proposed Partial Charges for Methanol (in fundamental charge units):")
print(f"Carbon: {charge_C:.4f}")
print(f"Methyl H1: {charge_H_methyl:.4f}")
print(f"Methyl H2: {charge_H_methyl:.4f}")
print(f"Methyl H3: {charge_H_methyl:.4f}")
print(f"Oxygen: {charge_O:.4f}")
print(f"Hydroxyl hydrogen: {charge_H_hydroxyl:.4f}")
print("\nVerifying that the molecule is charge-neutral:")
print(f"({charge_C:.4f}) + 3 * ({charge_H_methyl:.4f}) + ({charge_O:.4f}) + ({charge_H_hydroxyl:.4f}) = {total_charge:.4f}")

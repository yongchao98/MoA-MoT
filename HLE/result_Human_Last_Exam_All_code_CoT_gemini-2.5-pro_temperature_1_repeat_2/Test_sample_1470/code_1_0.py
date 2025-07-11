# Propose a reasonable set of partial charges for a new all-atom model of methanol (CH3OH).
# The charges are based on chemical principles and constraints of molecular simulations.

# Define the partial charges for each atom type in elementary charge units.
charge_C = 0.1450
charge_H_methyl = 0.0400
charge_O = -0.6830
charge_H_hydroxyl = 0.4180

# The methanol molecule has 1 Carbon, 3 methyl Hydrogens, 1 Oxygen, and 1 hydroxyl Hydrogen.
total_charge = charge_C + (3 * charge_H_methyl) + charge_O + charge_H_hydroxyl

# Print the proposed charge assignments for clarity.
print("Proposed Partial Charges for Methanol (CH3OH):")
print(f"Carbon: \t{charge_C:.4f}")
print(f"Methyl H1: \t{charge_H_methyl:.4f}")
print(f"Methyl H2: \t{charge_H_methyl:.4f}")
print(f"Methyl H3: \t{charge_H_methyl:.4f}")
print(f"Oxygen: \t{charge_O:.4f}")
print(f"Hydroxyl H: \t{charge_H_hydroxyl:.4f}")
print("-" * 35)

# Verify and print that the total charge of the molecule is zero.
print("Verification of Charge Neutrality:")
print(f"{charge_C:.4f} (C) + 3*({charge_H_methyl:.4f}) (H) + {charge_O:.4f} (O) + {charge_H_hydroxyl:.4f} (H) = {total_charge:.4f}")
import math

# Proposing a novel and reasonable set of partial charges for methanol (CH3OH)
# These charges are inspired by common force fields but are not a direct copy.

# Principles:
# 1. Oxygen is the most electronegative atom.
# 2. The hydroxyl hydrogen is highly positive due to bonding with oxygen.
# 3. Methyl hydrogens are slightly positive.
# 4. The whole molecule must be neutral (sum of charges = 0).
# 5. The three methyl hydrogens are equivalent.

# Charge assignments in elementary charge units
charge_C = 0.151
charge_H_methyl = 0.038
charge_O = -0.685
charge_H_hydroxyl = 0.420

# The molecule has 1 Carbon, 3 methyl Hydrogens, 1 Oxygen, and 1 hydroxyl Hydrogen.
net_charge = charge_C + (3 * charge_H_methyl) + charge_O + charge_H_hydroxyl

print("Proposed Partial Charges for Methanol (CH3OH):")
print("-" * 45)
print(f"Carbon (C):           {charge_C: .4f}")
print(f"Methyl Hydrogen (H):    {charge_H_methyl: .4f}")
print(f"Oxygen (O):           {charge_O: .4f}")
print(f"Hydroxyl Hydrogen (H):  {charge_H_hydroxyl: .4f}")
print("-" * 45)

# Verify that the net charge is zero. Floating point math can be imprecise,
# so we check if the sum is very close to zero.
print("\nVerification of Net Charge Neutrality:")
print(f"Net Charge = (C) + 3*(H_methyl) + (O) + (H_hydroxyl)")
print(f"Net Charge = {charge_C} + 3*({charge_H_methyl}) + ({charge_O}) + {charge_H_hydroxyl}")

# Using f-string for formatted output of the final sum
print(f"Sum = {net_charge:.4f}")

if math.isclose(net_charge, 0.0, abs_tol=1e-9):
    print("\nThe partial charges sum to zero. The molecule is neutral as required.")
else:
    print("\nWarning: The partial charges do not sum to zero!")

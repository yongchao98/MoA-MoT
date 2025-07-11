# Define the partial charges for each atom in the methanol molecule.
# These values are chosen based on chemical principles of electronegativity,
# molecular symmetry, and overall charge neutrality.

c_charge = 0.1450
h_methyl1_charge = 0.0400
h_methyl2_charge = 0.0400
h_methyl3_charge = 0.0400
o_charge = -0.6830
h_hydroxyl_charge = 0.4180

# Print the proposed charge assignments for the new model.
print("Proposed Partial Charges for Methanol (in fundamental charge units, e):")
print(f"Carbon:           {c_charge:.4f}")
print(f"Methyl H1:        {h_methyl1_charge:.4f}")
print(f"Methyl H2:        {h_methyl2_charge:.4f}")
print(f"Methyl H3:        {h_methyl3_charge:.4f}")
print(f"Oxygen:           {o_charge:.4f}")
print(f"Hydroxyl Hydrogen:  {h_hydroxyl_charge:.4f}")
print("-" * 30)

# Verify that the total charge of the molecule is zero.
total_charge = c_charge + h_methyl1_charge + h_methyl2_charge + h_methyl3_charge + o_charge + h_hydroxyl_charge

# Print the summation to explicitly show all numbers in the final equation.
print("Verification of Net Neutrality:")
print(f"Total charge = {c_charge:.4f} + {h_methyl1_charge:.4f} + {h_methyl2_charge:.4f} + {h_methyl3_charge:.4f} + ({o_charge:.4f}) + {h_hydroxyl_charge:.4f} = {total_charge:.4f}")

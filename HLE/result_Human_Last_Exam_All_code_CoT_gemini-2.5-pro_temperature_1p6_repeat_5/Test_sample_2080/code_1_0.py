# Plan: Calculate the second osmotic virial coefficient from steric-only behavior
# for a monoclonal antibody (mAb).

# Step 1: Define the standard partial specific volume for an mAb.
# This value is a physical constant for the protein and is given in mL/g.
# A typical value for mAbs is 0.73 mL/g.
partial_specific_volume = 0.73 # in mL/g

# Step 2: Use the hard-sphere model to calculate the steric-only contribution (B22,steric).
# The formula is B22,steric = 4 * partial_specific_volume.
# This term is always positive, representing hard-body repulsion.
factor = 4
b22_steric = factor * partial_specific_volume

# Step 3: Use the provided experimental total coefficient to provide context.
# B22,total = B22,steric + B22,electrostatic
b22_total = -7.585 # in mL/g
b22_electrostatic = b22_total - b22_steric

# Step 4: Print the calculation and the results in a clear format.
# The main request is for the steric-only coefficient.
print("--- Calculation of the Steric-Only Second Osmotic Virial Coefficient (B22,steric) ---")
print(f"Based on the hard-sphere model, B22,steric is calculated as 4 times the partial specific volume (v_sp).")
print(f"Using a standard v_sp for monoclonal antibodies of {partial_specific_volume} mL/g:")
print("\nFinal Equation:")
print(f"B22,steric = {factor} * {partial_specific_volume} mL/g")
print(f"B22,steric = {b22_steric:.3f} mL/g")


print("\n--- Context using the provided experimental data ---")
print(f"The given total coefficient (B22,total) is {b22_total} mL/g.")
print("This total value is the sum of steric and electrostatic contributions.")
print("B22,electrostatic = B22,total - B22,steric")
print(f"B22,electrostatic = {b22_total} - {b22_steric:.3f}")
print(f"B22,electrostatic = {b22_electrostatic:.3f} mL/g")
print("\nThe large negative electrostatic term confirms the presence of strong attractive interactions at pH 5 and 15 mM ionic strength.")

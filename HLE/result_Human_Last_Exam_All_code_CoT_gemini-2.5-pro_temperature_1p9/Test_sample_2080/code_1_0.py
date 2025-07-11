import sys
# Define the typical physical constants for a monoclonal antibody (mAb).
# The partial specific volume represents the volume occupied by the protein itself.
partial_specific_volume = 0.73  # Units: mL/g

# Hydration represents the volume of water tightly bound to the protein surface.
hydration = 0.4  # Units: mL/g

# Step 1: Calculate the hydrated specific volume of the mAb.
# This is the effective volume of the protein plus its hydration shell in solution.
hydrated_specific_volume = partial_specific_volume + hydration

# Step 2: Calculate the steric-only second osmotic virial coefficient (B22,steric).
# For a hard-sphere model, this coefficient (in units of mL/g) is four times
# the hydrated specific volume. This term represents the repulsive interactions
# due to the physical size of the molecules.
b22_steric = 4 * hydrated_specific_volume

# Print the explanation and the full calculation.
print("The steric-only contribution to the second osmotic virial coefficient (B22,steric) is calculated based on a hard-sphere model.")
print("The formula, for B22 in units of mL/g, is:")
print("B22,steric = 4 * (Partial Specific Volume + Hydration)")
print("\nUsing typical literature values for a monoclonal antibody:")
print(f"B22,steric = 4 * ({partial_specific_volume} mL/g + {hydration} mL/g)")
print(f"B22,steric = 4 * {hydrated_specific_volume:.2f} mL/g")
print(f"\nThe calculated second osmotic virial coefficient from steric-only behavior is: {b22_steric:.3f} mL/g")

# Use a special format to output the final numerical answer for automated reading.
# We redirect it to stderr to keep the explanation clean in stdout.
print(f"<<<{b22_steric:.3f}>>>", file=sys.stderr)
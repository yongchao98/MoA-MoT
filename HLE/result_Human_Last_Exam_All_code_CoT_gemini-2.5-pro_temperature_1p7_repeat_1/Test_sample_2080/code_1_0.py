# The task is to calculate the second osmotic virial coefficient from steric-only behavior for a monoclonal antibody.

# Define the typical partial specific volume for a monoclonal antibody in mL/g.
partial_specific_volume = 0.73

# The factor used to estimate the steric-only virial coefficient from the partial specific volume.
steric_factor = 8

# Calculate the steric-only second osmotic virial coefficient (B_22_steric).
b22_steric = steric_factor * partial_specific_volume

# Print an explanation and the final equation with the calculated result.
print("The steric-only contribution to the second osmotic virial coefficient (B_22_steric) is estimated from the partial specific volume (v_bar).")
print("The formula is: B_22_steric = 8 * v_bar")
print(f"Using a typical partial specific volume for an antibody of {partial_specific_volume} mL/g, the calculation is:")
print(f"B_22_steric = {steric_factor} * {partial_specific_volume} mL/g")
print(f"The result is: {b22_steric:.3f} mL/g")
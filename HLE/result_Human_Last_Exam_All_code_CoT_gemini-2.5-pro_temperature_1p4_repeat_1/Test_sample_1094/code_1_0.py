# This script presents the formula for the normalized AC loss in a superconducting
# elliptic bar carrying a transport current, as requested.

# The normalized AC loss is defined as L_norm = 2*pi*Q / (mu_0 * Ic^2),
# where Q is the loss per cycle per unit length, mu_0 is the magnetic permeability
# of free space, and Ic is the critical current of the bar.

# The variable 'i' represents the normalized current amplitude, i = Im / Ic,
# for the condition i < 1.

# The formula, derived from the critical state model (Norris model), is printed below.
# Note that 'ln' denotes the natural logarithm.
final_formula = "2 * ((1 - i) * ln(1 - i) + (1 + i) * ln(1 + i) - i^2)"

print("The normalized loss per cycle per unit length, 2*pi*Q/(mu_0*Ic^2), as a function of i = Im/Ic is:")
print(final_formula)
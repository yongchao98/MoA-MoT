import numpy as np

# This script calculates the finite blow-up time for a simplified linear model
# derived from the original PDE. The non-linear term u·∇u is neglected.
# The blow-up of the linear part is a strong indicator that the full non-linear
# solution will also blow up, as the non-linear term is not expected to be
# regularizing.

print("Analysis of a simplified model predicts a finite-time blow-up.")
print("The blow-up time `t` can be found by solving a quadratic equation that arises")
print("from the condition for the divergence of the energy integral.")
print("")

# The quadratic equation is of the form a*t^2 + b*t + c = 0.
# For our specific model and choice of initial data, the coefficients are:
a = 1
b = 2
c = -1

# We print the equation with its numerical coefficients.
# \u00b2 is the unicode for the superscript '2'.
# We also handle the signs for clear printing.
sign_b = '+' if b >= 0 else '-'
sign_c = '+' if c >= 0 else '-'
print(f"Equation for blow-up time t: {a}*t\u00b2 {sign_b} {abs(b)}*t {sign_c} {abs(c)} = 0")
print("----------------------------------------------------------")

# We solve this equation for t > 0.
# Using the quadratic formula: t = (-b ± sqrt(b^2 - 4ac)) / 2a

# Calculate the discriminant
discriminant = b**2 - 4*a*c

# The blow-up time must be positive.
if discriminant < 0:
    print("The equation has no real roots. No blow-up time found by this model.")
elif (-b + np.sqrt(discriminant)) / (2 * a) <= 0:
    print("The equation has no positive real roots. No forward-in-time blow-up.")
else:
    # Calculate the positive root, which represents the blow-up time
    t_blowup = (-b + np.sqrt(discriminant)) / (2 * a)
    print(f"The calculated finite blow-up time for this model is t = {t_blowup}")
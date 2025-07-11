# The problem is to find the normalized AC loss in an elliptical superconductor.
# The normalized loss is L(i) = 2*pi*Q/(mu_0 * Ic^2) as a function of i = Im/Ic.
# Based on the Norris formula for transport current loss in an ellipse for i < 1,
# the final equation is derived as L(i) = 2*(1-i)*ln(1-i) + 2*i - i^2.

# The final equation for normalized loss is of the form:
# L(i) = C1 * (1-i)*ln(1-i) + C2*i - C3*i^2

# Define the coefficients of this equation.
coeff_term1 = 2
coeff_term2 = 2
coeff_term3 = 1

# Print the final equation for the user, explicitly showing each coefficient.
print("The normalized loss per cycle, L(i) = 2*pi*Q/(mu_0*Ic^2), as a function of i=Im/Ic for i<1 is:")
print(f"L(i) = {coeff_term1} * (1-i)*ln(1-i) + {coeff_term2}*i - {coeff_term3}*i^2")
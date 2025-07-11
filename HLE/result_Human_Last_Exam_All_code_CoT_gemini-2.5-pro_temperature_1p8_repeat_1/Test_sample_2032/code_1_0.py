# All calculations are done using fractions to maintain precision.
from fractions import Fraction

# E[Y]
E_Y = Fraction(1, 2)

# E[Y^2]
E_Y_squared = Fraction(19, 60)

# Var(Y) = E[Y^2] - (E[Y])^2
variance_Y = E_Y_squared - E_Y**2

# Output the components of the variance calculation as requested.
# Using .limit_denominator() to show the fractions in their simplest form.
E_Y_squared_str = E_Y_squared.limit_denominator()
E_Y_str = E_Y.limit_denominator()
variance_Y_str = variance_Y.limit_denominator()

# Print the final equation
print(f"E[Y^2] = {E_Y_squared_str}")
print(f"E[Y] = {E_Y_str}")
print(f"Var(Y) = E[Y^2] - (E[Y])^2")
print(f"Var(Y) = {E_Y_squared_str} - ({E_Y_str})^2 = {E_Y_squared_str} - {E_Y**2} = {variance_Y_str}")

# For a numerical answer, we can convert the fraction to a float
numerical_variance = float(variance_Y)
print(f"The numerical value of the variance is: {numerical_variance}")
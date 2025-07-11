import math
from fractions import Fraction

# The problem is to determine the prefactor c_n for the fully f-connected
# diagram in the virial expansion for the coefficient B_n.

# From the analysis of the Mayer cluster expansion, the formula for c_n is
# derived to be:
# c_n = -(n - 1) / n!

# The code will demonstrate this calculation for the specific case of n = 4.

n = 4

# Calculate the components of the formula
# The numerator is -(n-1)
numerator_val = n - 1

# The denominator is n! (n factorial)
denominator_val = math.factorial(n)

# The prefactor c_n is the fraction of these two numbers.
# We use the Fraction class from the 'fractions' module to get an exact result.
c_n_fraction = Fraction(-numerator_val, denominator_val)

# As requested, we print each number involved in the final equation.
print(f"Calculation of the prefactor c_n for n = {n}:")
print(f"The general formula is c_n = - (n - 1) / n!")
print("\nPlugging in the numbers for n = 4:")
print(f"The term (n - 1) is: ({n} - 1) = {numerator_val}")
print(f"The term n! is: {n}! = {denominator_val}")
print("\nThus, the final prefactor is:")
print(f"c_{n} = -{numerator_val} / {denominator_val} = {c_n_fraction}")
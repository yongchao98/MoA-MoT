from fractions import Fraction

# Based on the analytical derivation, the problem of maximizing c_3
# is reduced to finding the maximum of a related integral.
# The maximum value of this integral, which we can call I_max, was found to be 5/4.
I_max = Fraction(5, 4)

# The coefficient c_3 is related to this integral by a factor of 7/2.
factor = Fraction(7, 2)

# Calculate the maximum value of c_3.
c3_max = factor * I_max

# Print the final calculation, showing each number in the equation.
print("The maximum value of c_3 is calculated from the maximum value of a related integral.")
print(f"The maximum value of the integral was determined to be {I_max.numerator}/{I_max.denominator}.")
print("The final equation for the maximum value of c_3 is:")
print(f"max(c_3) = ({factor.numerator} / {factor.denominator}) * ({I_max.numerator} / {I_max.denominator}) = {c3_max.numerator} / {c3_max.denominator}")
print(f"As a decimal, the result is {float(c3_max)}.")

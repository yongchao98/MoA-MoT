# In the context of a G4-theoretical framework, the upper critical dimension is 4.
# For dimensions d >= 4, mean-field theory is exact. According to mean-field
# theory, the critical exponent nu (ν) for the correlation length is precisely 1/2.
#
# This script calculates this value and displays the equation.

# Define the numerator and denominator from the mean-field theory result.
numerator = 1
denominator = 2

# Calculate the precise value of the critical exponent ν.
nu_value = numerator / denominator

# As requested, print each number in the final equation.
# The equation shows the derivation of ν from its mean-field value.
print(f"ν = {numerator} / {denominator} = {nu_value}")
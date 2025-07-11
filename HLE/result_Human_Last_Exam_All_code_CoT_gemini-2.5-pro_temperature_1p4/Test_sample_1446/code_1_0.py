# The problem asks for the critical exponent nu (ν) in a G₄-theoretical
# framework. This is interpreted as the standard φ⁴ (phi-four) model of
# critical phenomena.
#
# The value of ν depends on the spatial dimension d. The request for a
# single "precise value" points to the mean-field theory result, which is
# exact at and above the upper critical dimension (d_c = 4).
#
# In this mean-field regime, the critical exponent ν is exactly 1/2.

# The equation for nu is ν = 1 / 2.
# We define the numerator and denominator to show each number.
numerator = 1
denominator = 2

# Calculate the precise value of the exponent
nu_value = numerator / denominator

# Print the components of the equation and the final result.
print(f"The equation for the mean-field critical exponent is: nu = {numerator} / {denominator}")
print(f"The precise value of the critical exponent nu (ν) is: {nu_value}")
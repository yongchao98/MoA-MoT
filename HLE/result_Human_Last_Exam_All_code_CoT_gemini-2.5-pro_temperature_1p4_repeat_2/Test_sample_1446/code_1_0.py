# This script calculates the value of the critical exponent nu (ν)
# for a G₄ (phi-4) theory at or above its upper critical dimension (d >= 4).
# In this regime, the system's behavior is described by mean-field theory,
# which provides an exact value for the exponent.

# In mean-field theory, the exponent ν is given by the fraction 1/2.
numerator = 1
denominator = 2

# Calculate the precise decimal value.
nu_value = numerator / denominator

# Print the result in an equation format as requested, showing each number.
print(f"The critical exponent ν is calculated as:")
print(f"{numerator} / {denominator} = {nu_value}")
